/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For RHS update
#include "BoostedIsotropicKerr.hpp"
//#include "FixedBGEvolution.hpp"
#include "FAKEMatterCCZ4RHS.hpp"
#include "FourthOrderDerivatives.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComplexPotential.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "ComplexScalarField.hpp"
#include "FixedBGEnergyAndMomFlux.hpp"
#include "FixedBGMomAndSource.hpp"
#include "FluxExtraction.hpp"
#include "ScalarConstant.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for scalar field
    SetValue set_zero(0.0);
    BoostedIsotropicKerr boosted_bh(m_p.bg_params, m_dx);
    ScalarConstant initial_sf(m_p.scalar_amplitude, m_p.scalar_mass, m_p.center,
                              m_p.bg_params, m_dx);
    auto compute_pack = make_compute_pack(set_zero, boosted_bh);


    BoxLoops::loop(compute_pack, m_state_new, m_state_new,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential>(
            m_dx, m_p.center, m_p.bg_params),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    int min_level = 0;
    bool calculate_quantities = at_level_timestep_multiple(min_level);

    if (calculate_quantities)
    {
        fillAllGhosts();
        ComplexPotential potential(m_p.scalar_mass);
        ScalarFieldWithPotential scalar_field(potential);
        BoostedIsotropicKerr boosted_bh(m_p.bg_params, m_dx);
        FixedBGMomAndSource<ScalarFieldWithPotential, BoostedIsotropicKerr>
            densities(scalar_field, boosted_bh, m_dx, m_p.center);
        FixedBGEnergyAndMomFlux<ScalarFieldWithPotential, BoostedIsotropicKerr>
            fluxes(scalar_field, boosted_bh, m_dx, m_p.center);
        BoxLoops::loop(make_compute_pack(densities, fluxes), m_state_new,
                       m_state_diagnostics, SKIP_GHOST_CELLS);

        // excise within horizon, no simd
        
        BoxLoops::loop(
            ExcisionDiagnostics<ScalarFieldWithPotential, BoostedIsotropicKerr>(
                m_dx, m_p.center, boosted_bh, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
        
    }

    // write out the integral after each coarse timestep
    if (m_level == 0)
    {
        bool first_step = (m_time == m_dt);
        // integrate the densities and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double xSource_sum = amr_reductions.sum(c_xSource);
        double ySource_sum = amr_reductions.sum(c_ySource);
        double zSource_sum = amr_reductions.sum(c_zSource);
        double xMom_sum = amr_reductions.sum(c_xMom);
        double yMom_sum = amr_reductions.sum(c_yMom);
        double zMom_sum = amr_reductions.sum(c_zMom);
        double rho_sum = amr_reductions.sum(c_rho);

        SmallDataIO integral_file("SourceXMomRhoInts", m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();

        std::vector<double> data_for_writing = {xSource_sum, ySource_sum, zSource_sum,
                                 xMom_sum, yMom_sum, zMom_sum, rho_sum};

        // write data
        if (first_step)
        {
            integral_file.write_header_line({"xSource", "ySource", "zSource",
                                             "x-Mom", "y-Mom", "z-Mom", "rho"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // Now refresh the interpolator and do the interpolation
        bool fill_ghosts = false;
        m_gr_amr.m_interpolator->refresh(fill_ghosts);
        m_gr_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                        Interval(c_xMdot, c_Edot));
        FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                     m_restart_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }
}

// Things to do before a plot level - need to calculate the Stress
void ScalarFieldLevel::prePlotLevel() {}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate right hand side with matter_t = ScalarField
    // and background_t = BoostedBH
    // RHS for non evolution vars is zero, to prevent undefined values
    ComplexPotential potential(m_p.scalar_mass);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedIsotropicKerr boosted_bh(m_p.bg_params, m_dx);
    MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    /*
    FixedBGEvolution<ScalarFieldWithPotential, BoostedIsotropicKerrFixedBG> my_evolution(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);

    BoxLoops::loop(my_evolution, a_soln, a_rhs, SKIP_GHOST_CELLS);
    */

    // Do excision within horizon
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential>(
            m_dx, m_p.center, m_p.bg_params),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());

    
}

// Note that for the fixed grids this only happens on the initial timestep
// simd is disabled to allow simpler use of logical operators
void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    //pout() << "m_P center: " << m_p.center[0] << m_p.center[1] << m_p.center[2] << endl;
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion, disable_simd());
}

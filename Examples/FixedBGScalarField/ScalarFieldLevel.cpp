/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"

// For RHS update
#include "BoostedBHFixedBG.hpp"
#include "Excision.hpp"
#include "FixedBGEvolution.hpp"

// For density calculation
#include "FixedBGDensity.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "FixedBGScalarField.hpp"
#include "Potential.hpp"
#include "ScalarGaussian.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for scalar field
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    ScalarGaussian initial_sf(m_p.initial_params, m_dx);
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boosted_bh, initial_sf),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
    // Calculate matter density function
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    BoxLoops::loop(FixedBGDensity<ScalarFieldWithPotential, BoostedBHFixedBG>(
                       scalar_field, boosted_bh, m_dx, m_p.center),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate right hand side with matter_t = ScalarField
    // and background_t = BoostedBH
    // RHS for non evolution vars is zero, to prevent undefined values
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ScalarFieldWithPotential, BoostedBHFixedBG> my_evolution(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);
    SetValue set_static_rhs_zero(0.0, Interval(c_chi, c_rho));
    auto compute_pack = make_compute_pack(my_evolution, set_static_rhs_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

    // Do excision within horizon
    BoxLoops::loop(Excision<ScalarFieldWithPotential, BoostedBHFixedBG>(
                       m_dx, m_p.center, boosted_bh),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Specify if you want any plot files to be written, with which vars
void ScalarFieldLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = {c_phi, c_chi};
}

// Note that for the fixed grids this only happens on the initial timestep
// simd is disabled to allow simpler use of logical operators
void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion, disable_simd());
}

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

// For momentum flux calculation
#include "FixedBGStress.hpp"

// For tag cells
#include "ExtractionFixedGridsTaggingCriterion.hpp"
//#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "FixedBGComplexScalarField.hpp"
#include "ComplexPotential.hpp"
#include "ScalarConstant.hpp"
#include "SetValue.hpp"

//#include "WeylExtraction2.hpp"
#include "StressExtraction.hpp"
#include "StressExtraction_halves.hpp"
#include "CustomExtraction.hpp"

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
    ScalarConstant initial_sf(m_p.scalar_amplitude, m_p.scalar_mass, m_p.center,
			      m_p.bg_params, m_dx);
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boosted_bh, initial_sf),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    
    //setup the output file
    SmallDataIO integral_file(m_p.integral_filename, m_dt, m_time,
                              m_restart_time, SmallDataIO::APPEND, true);
    std::vector<std::string> header_strings = {"Source","Xmom"};
    integral_file.write_header_line(header_strings);
}

void ScalarFieldLevel::specificPostTimeStep()
{
  CH_TIME("ScalarFieldLevel::specificPostTimeStep");
  if (m_p.activate_extraction == 1)
    {
      // Populate the Stress values on the grid
      fillAllGhosts();
      ComplexPotential potential(m_p.scalar_mass);
      ScalarFieldWithPotential scalar_field(potential);
      BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
      BoxLoops::loop(FixedBGStress<ScalarFieldWithPotential, 
		     BoostedBHFixedBG>(scalar_field, boosted_bh, m_dx, m_p.center),
                     m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
      
      // write out the integral after each coarse timestep
      if (m_level == 0)
	{
	  // integrate the densities and write to a file
	  double Source_sum = m_gr_amr.compute_sum(c_Source, m_p.coarsest_dx);
	  double Xmom_sum = m_gr_amr.compute_sum(c_Xmom, m_p.coarsest_dx);
	  
	  SmallDataIO integral_file(m_p.integral_filename, m_dt, m_time,
				    m_restart_time, SmallDataIO::APPEND, false);
	  // remove any duplicate data if this is post restart
	  integral_file.remove_duplicate_time_data();
	  std::vector<double> data_for_writing = {Source_sum, Xmom_sum};
	  // write data
	  integral_file.write_time_data_line(data_for_writing);
	  //	}
	  //if (m_level == m_p.extraction_params.min_extraction_level)
	  //{	  
	  // Now refresh the interpolator and do the interpolation
	  m_gr_amr.m_interpolator->refresh();
	  StressExtraction my_extraction(m_p.extraction_params, m_dt, m_time, m_restart_time);
	  my_extraction.execute_query(m_gr_amr.m_interpolator);

          m_gr_amr.m_interpolator->refresh();
          StressExtraction_halves my_extraction_halves(m_p.extraction_params, m_dt, m_time, m_restart_time);
          my_extraction_halves.execute_query(m_gr_amr.m_interpolator);
	  //} 
	  //if (m_level == 2)
	  //{
	  m_gr_amr.m_interpolator->refresh();
          int num_points = 200;
          CustomExtraction my_extraction_phi(c_phi_Re, num_points, m_p.L, m_p.center,
					     m_dt, m_time);
          my_extraction_phi.execute_query(m_gr_amr.m_interpolator,
                                          m_p.extraction_filename);

          m_gr_amr.m_interpolator->refresh();
          CustomExtraction my_extraction_rho(c_rho, num_points, m_p.L, m_p.center,
                                             m_dt, m_time);
          my_extraction_rho.execute_query(m_gr_amr.m_interpolator,
                                          m_p.extraction_filename2);
	}
    }
}
// Things to do before a plot level - need to calculate the Stress
void ScalarFieldLevel::prePlotLevel()
{
  //  fillAllGhosts();
  //if (m_p.activate_extraction == 1)
  //  {
  //    ComplexPotential potential(m_p.scalar_mass);
  //    ScalarFieldWithPotential scalar_field(potential);
  //    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
  //    BoxLoops::loop(FixedBGStress<ScalarFieldWithPotential, BoostedBHFixedBG>(
  //	                 scalar_field, boosted_bh, m_dx, m_p.center),
  //		     m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
  //}
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate right hand side with matter_t = ScalarField
    // and background_t = BoostedBH
    // RHS for non evolution vars is zero, to prevent undefined values
    ComplexPotential potential(m_p.scalar_mass);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ScalarFieldWithPotential, BoostedBHFixedBG> my_evolution(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);
    SetValue set_static_rhs_zero(0.0, Interval(c_chi,c_dArea));
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
  plot_states = {c_phi_Re, c_phi_Im, c_rho, c_Source, c_Xmom, c_Stress, c_dArea};
}

// Note that for the fixed grids this only happens on the initial timestep
// simd is disabled to allow simpler use of logical operators
void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
  BoxLoops::loop(ExtractionFixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center,  m_p.extraction_params),
		 current_state, tagging_criterion, disable_simd());
  //BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
  //		 current_state, tagging_criterion, disable_simd());
}

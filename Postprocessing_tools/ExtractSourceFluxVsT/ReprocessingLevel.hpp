/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef REPROCESSINGLEVEL_HPP_
#define REPROCESSINGLEVEL_HPP_

#include "GRAMRLevel.hpp"
#include "ExcisionSource.hpp"

#include "BoxLoops.hpp"
#include "NanCheck.hpp"

// For momentum flux calculation
#include "FixedBGStress.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "FixedBGComplexScalarField.hpp"
#include "ComplexPotential.hpp"
#include "ScalarConstant.hpp"
#include "SetValue.hpp"

#include "StressExtraction.hpp"

class ReprocessingLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ReprocessingLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    typedef FixedBGComplexScalarField<ComplexPotential> ScalarFieldWithPotential;

    // initialize data
    virtual void initialData() { m_state_new.setVal(0.); }

    void postRestart()
    {
      if (m_p.activate_extraction == 1)
        {
          // Populate the Stress values on the grid
          fillAllGhosts();
          ComplexPotential potential(m_p.scalar_mass);
          ScalarFieldWithPotential scalar_field(potential);
          BoostedIsotropicBHFixedBG boosted_bh(m_p.bg_params, m_dx);
	  //BoxLoops::loop(FixedBGStress<ScalarFieldWithPotential,
	  //		 BoostedIsotropicBHFixedBG>(scalar_field, boosted_bh, m_dx, m_p.center),
	  //		 m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
	  BoxLoops::loop(ExcisionSource<ScalarFieldWithPotential, BoostedIsotropicBHFixedBG>(
        		 m_dx, m_p.center, boosted_bh, true),
			 m_state_new, m_state_new, SKIP_GHOST_CELLS,
			 disable_simd());

          if (m_level == 9)
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

              // Now refresh the interpolator and do the interpolation
              //m_gr_amr.m_interpolator->refresh();
	      AMRInterpolator<Lagrange<4>> my_interpolator(m_gr_amr, m_p.origin, m_p.dx, m_p.verbosity);
	      my_interpolator.refresh();
              StressExtraction my_extraction(m_p.extraction_params, m_dt, m_time, m_restart_time);
              //my_extraction.execute_query(m_gr_amr.m_interpolator);
	      my_extraction.execute_query(&my_interpolator); 
            }
        }
        pout() << "The time is " << m_time << " on level " << m_level 
               << ". Your wish is my command." << endl;
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* REPROCESSINGLEVEL_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "CCZ4.hpp"
#include "BoundaryConditions.hpp"

// Problem specific includes:
#include "BoostedBHFixedBG.hpp"
#include "Potential.hpp"
#include "ScalarConstant.hpp"

struct extraction_params_t
{
  int num_extraction_radii;
  std::vector<double> extraction_radii;
  std::array<double, CH_SPACEDIM> extraction_center;
  int num_points_phi;
  int num_points_theta;
  int num_modes;
  std::vector<std::pair<int, int>> modes; // l = first, m = second                                                                                                         
  std::vector<int> extraction_levels;
  bool write_extraction;
  int min_extraction_level;
};

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for nan check
        pp.load("nan_check", nan_check, 1);

        // Initial and SF data
        pp.load("scalar_mass", potential_params.scalar_mass);
        pp.load("scalar_amplitude", initial_params.amplitude);
        //pp.load("scalar_width", initial_params.width);
        //pp.load("scalar_center", initial_params.center, center);
        pp.load("sigma", sigma);

        // Background boosted bh data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);

        pp.load("activate_extraction", activate_extraction, 0);

	// extraction params
        dx.fill(coarsest_dx);
        origin.fill(coarsest_dx / 2.0);

        // Extraction params
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii,
                1);
        // Check for multiple extraction radii, otherwise load single
        // radius/level (for backwards compatibility).
        if (pp.contains("extraction_levels"))
	  {
            pp.load("extraction_levels", extraction_params.extraction_levels,
                    extraction_params.num_extraction_radii);
	  }
        else
	  {
            pp.load("extraction_level", extraction_params.extraction_levels, 1,
                    0);
	  }
        if (pp.contains("extraction_radii"))
	  {
            pp.load("extraction_radii", extraction_params.extraction_radii,
                    extraction_params.num_extraction_radii);
	  }
        else
	  {
            pp.load("extraction_radius", extraction_params.extraction_radii, 1,
                    0.1);
	  }
        pp.load("num_points_phi", extraction_params.num_points_phi, 2);
        pp.load("num_points_theta", extraction_params.num_points_theta, 4);
        pp.load("extraction_center", extraction_params.extraction_center,
                center);
        if (pp.contains("modes"))
	  {
            pp.load("num_modes", extraction_params.num_modes);
	    std::vector<int> extraction_modes_vect(2 *
                                                   extraction_params.num_modes);
            pp.load("modes", extraction_modes_vect,
                    2 * extraction_params.num_modes);
            extraction_params.modes.resize(extraction_params.num_modes);
            for (int i = 0; i < extraction_params.num_modes; ++i)
	      {
                extraction_params.modes[i].first = extraction_modes_vect[2 * i];
                extraction_params.modes[i].second =
		  extraction_modes_vect[2 * i + 1];
	      }
	  }
          else
	  {
            // by default extraction (l,m) = (2,0), (2,1) and (2,2)
            extraction_params.num_modes = 3;
            extraction_params.modes.resize(3);
            for (int i = 0; i < 3; ++i)
	      {
                extraction_params.modes[i].first = 2;
                extraction_params.modes[i].second = i;
	      }
	  }

        pp.load("write_extraction", extraction_params.write_extraction, false);

        // Work out the minimum extraction level
        auto min_extraction_level_it =
	  std::min_element(extraction_params.extraction_levels.begin(),
			   extraction_params.extraction_levels.end());
        extraction_params.min_extraction_level = *(min_extraction_level_it);
    }

    // Problem specific parameters
    double sigma;
    int nan_check;

    // Initial data for matter, metric and potential
    int activate_extraction;
    BoostedBHFixedBG::params_t bg_params;
    ScalarConstant::params_t initial_params;
    Potential::params_t potential_params;

    std::array<double, CH_SPACEDIM> origin,
          dx; // location of coarsest origin and dx

    // Collection of parameters necessary for the CCZ4 RHS and extraction 
    CCZ4::params_t ccz4_params;
    extraction_params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

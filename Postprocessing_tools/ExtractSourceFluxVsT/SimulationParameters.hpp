/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "FixedBGSimulationParametersBase.hpp"

// Problem specific includes:
#include "BoostedIsotropicBHFixedBG.hpp"
#include "ComplexPotential.hpp"
#include "ScalarConstant.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
public:
  SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
  {
    readParams(pp);
  }

    void readParams(GRParmParse &pp)
    {
        pp.get("num_files", num_files);
        pp.get("start_file", start_file);
        pp.get("checkpoint_interval", checkpoint_interval);

        // Initial and SF data
        pp.load("scalar_mass", scalar_mass);
        pp.load("scalar_amplitude", scalar_amplitude);
	//pp.load("scalar_center", initial_params.center, center);
	pp.load("integral_filename", integral_filename);
	//	pp.load("extraction_filename", extraction_filename);
	//pp.load("extraction_filename2", extraction_filename2);

        // Background boosted bh data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);

        pp.load("activate_extraction", activate_extraction, 0);
    }

    // Initial data for matter, metric and potential
    int activate_extraction;
    int num_files, start_file, checkpoint_interval;
    double scalar_mass, scalar_amplitude;
  std::string integral_filename; //extraction_filename, extraction_filename2;
      BoostedIsotropicBHFixedBG::params_t bg_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

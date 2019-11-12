/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "BoostedBHFixedBG.hpp"
#include "Potential.hpp"
#include "ScalarGaussian.hpp"

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
        pp.load("scalar_width", initial_params.width);
        pp.load("scalar_center", initial_params.center, center);
        pp.load("sigma", sigma);

        // Background boosted bh data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);
    }

    // Problem specific parameters
    double sigma;
    int nan_check;

    // Initial data for matter, metric and potential
    BoostedBHFixedBG::params_t bg_params;
    ScalarGaussian::params_t initial_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

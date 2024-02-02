/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "FixedBGSimulationParametersBase.hpp"
#include "SimulationParametersBase.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "BoostedIsotropicKerr.hpp"
#include "ComplexPotential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_length", regrid_length, L);

        // BH data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);
        pp.load("bh_spin", bg_params.spin);

        // Initial SF
        pp.load("scalar_amplitude", scalar_amplitude);
        pp.load("scalar_mass", scalar_mass);

        // Volume extraction radii
        pp.load("num_integral_r", num_integral_r, 1);
        pp.load("inner_r", inner_r, num_integral_r);
        pp.load("outer_r", outer_r, num_integral_r);
        pout() << "outer_r 0: " << outer_r[0] << endl;
        pout() << "outer_r 1: " << outer_r[1] << endl;
        pout() << "inner_r 0: " << inner_r[0] << endl;
        pout() << "inner_r 1: " << inner_r[1] << endl;
    }

    // Problem specific parameters
    double scalar_amplitude, scalar_mass, regrid_length;
    // Collection of parameters necessary for the sims
    BoostedIsotropicKerr::params_t bg_params;
    double G_Newton;
    int num_integral_r;
    std::vector<double> inner_r, outer_r;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

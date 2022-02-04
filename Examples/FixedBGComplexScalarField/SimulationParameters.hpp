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
#include "ComplexPotential.hpp"
#include "SphericalExtraction.hpp"

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
        // for regridding
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);
        pp.load("regrid_length", regrid_length, L);

	// BH data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);

        // Initial SF
        pp.load("scalar_amplitude", scalar_amplitude);
        pp.load("scalar_mass", scalar_mass);

	// Volume extraction radii
        pp.load("inner_r", inner_r, 5.0);
        pp.load("outer_r", outer_r, 100.0 / scalar_mass);

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
        pp.load("num_points_theta", extraction_params.num_points_theta, 5);
        if (extraction_params.num_points_theta % 2 == 0)
        {
            extraction_params.num_points_theta += 1;
            pout() << "Parameter: num_points_theta incompatible with Simpson's "
                   << "rule so increased by 1.\n";
        }
        pp.load("extraction_center", extraction_params.center, center);
        pp.load("write_extraction", extraction_params.write_extraction, false);
    }

    // Problem specific parameters
    double  sigma, excision_width, regrid_length;
    double scalar_amplitude, scalar_mass;
    int nan_check;
    double inner_r, outer_r;
    //    std::array<double, CH_SPACEDIM> origin,
    //        dx; // location of coarsest origin and dx
    // Collection of parameters necessary for the sims
    BoostedBHFixedBG::params_t bg_params;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

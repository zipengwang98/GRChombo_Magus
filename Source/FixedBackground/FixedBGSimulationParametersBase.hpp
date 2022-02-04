/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGSIMULATIONPARAMETERSBASE_HPP_
#define FIXEDBGSIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "BoundaryConditions.hpp"
#include "SphericalExtraction.hpp"

using extraction_params_t = SphericalExtraction::params_t;

class FixedBGSimulationParametersBase : public ChomboParameters
{
  public:
    FixedBGSimulationParametersBase(GRParmParse &pp) : ChomboParameters(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

  private:
    void readParams(GRParmParse &pp)
    {
        // for nan check
        pp.load("nan_check", nan_check, 1);

        //Dissipation
	pp.load("sigma", sigma, 0.1);

        // Extraction params
        pp.load("activate_extraction", activate_extraction, false);

        if (activate_extraction)
	  {
            pp.load("num_extraction_radii",
                    extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("extraction_levels"))
	      {
                pp.load("extraction_levels",
                        extraction_params.extraction_levels,
                        extraction_params.num_extraction_radii);
	      }
            else
	      {
                pp.load("extraction_level", extraction_params.extraction_levels,
                        1, 0);
	      }

	    if (pp.contains("extraction_radii"))
	      {
		pp.load("extraction_radii", extraction_params.extraction_radii,
			extraction_params.num_extraction_radii);
	      }
	    else
	      {
		pp.load("extraction_radius", extraction_params.extraction_radii,
			1, 0.1);
	      }

	    pp.load("num_points_phi", extraction_params.num_points_phi, 2);
	    pp.load("num_points_theta", extraction_params.num_points_theta, 5);
	    if (extraction_params.num_points_theta % 2 == 0)
	      {
		extraction_params.num_points_theta += 1;
		pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
		       << "rule so increased by 1.\n";
	      }
	    pp.load("extraction_center", extraction_params.center, center);

	    if (pp.contains("modes"))
	      {
		pp.load("num_modes", extraction_params.num_modes);
		std::vector<int> extraction_modes_vect(
		   2 * extraction_params.num_modes);
		pp.load("modes", extraction_modes_vect,
			2 * extraction_params.num_modes);
		extraction_params.modes.resize(extraction_params.num_modes);
		for (int i = 0; i < extraction_params.num_modes; ++i)
		  {
                    extraction_params.modes[i].first =
		      extraction_modes_vect[2 * i];
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

	    pp.load("write_extraction", extraction_params.write_extraction,
		false);

	    std::string extraction_path;
	    pp.load("extraction_subpath", extraction_path, data_path);
	    if (!extraction_path.empty() && extraction_path.back() != '/')
	      extraction_path += "/";
	    if (output_path != "./" && !output_path.empty())
	      extraction_path = output_path + extraction_path;

	    extraction_params.data_path = data_path;
	    extraction_params.extraction_path = extraction_path;

	    // default names to Weyl extraction
	    pp.load("extraction_file_prefix",
		extraction_params.extraction_file_prefix,
		    std::string("Weyl4_extraction_"));
	    pp.load("integral_file_prefix",
		    extraction_params.integral_file_prefix,
		    std::string("Weyl4_mode_"));
	  }
    }
}

public:

    double sigma;
    int nan_check;

    // Collection of parameters necessary for the CCZ4 RHS and extraction 
    bool activate_extraction;
    extraction_params_t extraction_params;
};

#endif /* FIXEDBGSIMULATIONPARAMETERSBASE_HPP_ */

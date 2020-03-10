/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXTRACTIONFIXEDGRIDSTAGGINGCRITERION_HPP_
#define EXTRACTIONFIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "FixedBGSimulationParametersBase.hpp"
#include "Tensor.hpp"

class ExtractionFixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const extraction_params_t m_params;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    ExtractionFixedGridsTaggingCriterion(const double dx, const int a_level,
				       const double a_L, const std::array<double, CH_SPACEDIM> a_center,
                                       const extraction_params_t a_params)
      : m_dx(dx), m_deriv(dx),
        m_level(a_level), m_L(a_L), m_center(a_center),
        m_params(a_params){};

  void compute(Cell<double> current_cell) const
  {
    double criterion = 0.0;

    for (int iradius = 0; iradius < m_params.num_extraction_radii;
	   ++iradius)
	{
	  // regrid if within extraction level and not at required refinement                                                                                              
	  if (m_level < m_params.extraction_levels[iradius])
	    {
	      const Coordinates<double> coords_ext(current_cell, m_dx,
						   m_params.extraction_center);
	      const double r = coords_ext.get_radius();
	      // add a 20% buffer to extraction zone so not too near to                                                                                                    
	      // boundary                                                                                                                                                  
	      bool regrid = r < 1.2 * m_params.extraction_radii[iradius];
	      if (regrid)
	      {
	       criterion = 100.0;
	      }
	    }
	}

    // make sure the inner part is regridded around the horizon
    for(int i_level = 0 ; i_level < 10; i_level ++)
      {
	if (m_level == i_level)
	  {
	    // take L as the length of full grid, so tag inner 1/2
	    // of it, which means inner \pm L/4
	    double ratio = pow(2.0, -(i_level + 2.0));
	    const Coordinates<double> coords_cen(current_cell, m_dx, m_center);
	    if (abs(coords_cen.x) < m_L*ratio &&
		abs(coords_cen.y) < m_L*ratio &&
		abs(coords_cen.z) < m_L*ratio)
	      {
		criterion = 100;
	      }
	  }
      }

    // Write back into the flattened Chombo box
    current_cell.store_vars(criterion, 0);
  }
};

#endif /* EXTRACTIONFIXEDGRIDSTAGGINGCRITERION_HPP_ */

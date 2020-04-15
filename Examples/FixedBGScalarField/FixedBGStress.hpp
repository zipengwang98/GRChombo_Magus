/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGSTRESS_HPP_
#define FIXEDBGSTRESS_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "CoordinateTransformations.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the grid
template <class matter_t, class background_t> class FixedBGStress
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

  public:
    FixedBGStress(matter_t a_matter, background_t a_background, double a_dx,
                   std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // get the metric vars from the background
        MetricVars<data_t> metric_vars;
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        m_background.compute_metric_background(metric_vars, coords);

        using namespace TensorAlgebra;
	using namespace CoordinateTransformations;
	//	const auto gamma = metric_vars.gamma;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
	const auto lapse = metric_vars.lapse;
	const auto shift = metric_vars.shift;

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t R = coords.get_radius();
	const auto gamma_spher = cartesian_to_spherical_LL(metric_vars.gamma, x, y, z);

	Tensor<1, data_t> sCart;
	sCart[0] = x/R;
	sCart[1] = y/R;
	sCart[2] = z/R;
	
	data_t sCart_norm = 0.0;
	FOR2(j, k)
	  {
	    sCart_norm += sCart[j]*sCart[k]*metric_vars.gamma[j][k];
	  }
	sCart_norm = sqrt(sCart_norm);
	Tensor<1, data_t> nCart;
	FOR1(i)
	{
	  nCart[i] = sCart[i]/sCart_norm;
	}

	Tensor<1, data_t> TUiDx;
        data_t Stress = 0.0;
	FOR1(i)
	{
	  TUiDx[i] = shift[i]*emtensor.Si[0]*nCart[i]/lapse;
	  FOR1(j)
	  {
	    TUiDx[i] += (metric_vars.gamma[i][j] - shift[i]*shift[j]/pow(lapse,2))*emtensor.Sij[j][0];
	  }
	}

	FOR1(i)
	{
	  Stress += TUiDx[i]*nCart[i];
	}
	Tensor<1, data_t> sSpher;
        sSpher[0] = 1.0;
        sSpher[1] = 0.0;
        sSpher[2] = 0.0;
	
        data_t sSpher_norm = sqrt(gamma_spher[0][0]);
        Tensor<1, data_t> nSpher; 
	FOR1(i)
        {
          nSpher[i] = sSpher[i]/sSpher_norm;
        }

	Tensor<2, data_t> Sigma;
	FOR2(i, j)
	  {
	    Sigma[i][j] = gamma_spher[i][j]; //+ gamma_spher[i][0]*gamma_spher[0][j];
	    FOR2(m,n)
	    {
	      Sigma[i][j] += gamma_spher[i][m] * nSpher[m] * gamma_spher[j][n] * nSpher[n];
	    }
	  }

	const data_t detSigma = Sigma[1][1] * Sigma[2][2] - Sigma[1][2] * Sigma[2][1];
	const data_t dArea = sqrt(Sigma[1][1] * Sigma[2][2] - Sigma[1][2] * Sigma[2][1]);
	//pout()<< "S10" << emtensor.Sij[0][1]<<emtensor.Sij[2][1]<<emtensor.Sij[2][3] <<endl;
	//pout()<< "S01" << emtensor.Sij[1][0]<<emtensor.Sij[1][2]<<emtensor.Sij[3][2] <<endl;
	//pout()<< "detgamma" << gamma_spher[1][1]*gamma_spher[2][2] - gamma_spher[1][2]*gamma_spher[2][1];
        // assign values of Momentum flux in output box
        current_cell.store_vars(Stress, c_Stress);
	current_cell.store_vars(dArea, c_dArea);
    }
};

#endif /* FIXEDBGSTRESS_HPP_ */

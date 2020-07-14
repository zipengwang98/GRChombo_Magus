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
	const data_t det_gamma = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t R = coords.get_radius();
	const auto gamma_spher = cartesian_to_spherical_LL(metric_vars.gamma, x, y, z);

	Tensor<1, data_t> si;
	si[0] = x/R;
	si[1] = y/R;
	si[2] = z/R;
	
	data_t si_norm = 0.0;
	FOR2(j, k)
	  {
	    si_norm += si[j]*si[k]*metric_vars.gamma[j][k];
	  }

	FOR1(i)
	{
	  si[i] = si[i]/sqrt(si_norm);
	}

	data_t Stress = 0.0;
	FOR1(i)
	{
          Stress += - metric_vars.lapse * si[i]*emtensor.Sij[i][0];
          FOR1(j)
          {
            Stress += metric_vars.gamma[i][j] * metric_vars.shift[j]
              * si[i] * emtensor.Si[0];
          }
	}
	
	const auto dArea = area_element_sphere(gamma_spher, x, y, z);
	data_t rho2 = simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
	data_t r2sintheta = sqrt(rho2) * R;
	
	data_t Stress2 = Stress * dArea / r2sintheta; 
	//pout()<< "si[0]" << si[0] << endl;
	//pout()<< "si[1]" << si[1] << endl;
	//pout()<< "si[2]" << si[2] << endl;
	//pout()<< "Sij[0][0]" << emtensor.Sij[0][0] <<endl;
	//pout()<< "Sij[1][0]" << emtensor.Sij[1][0] <<endl;
	//pout()<< "Sij[0][1]" << emtensor.Sij[0][1] <<endl;
	//pout()<< "Sij[2][0]" << emtensor.Sij[2][0] <<endl;
	//pout()<< "Sij[0][2]" << emtensor.Sij[0][2] <<endl;
	//pout()<< "Stress" << Stress << endl;
	//pout()<< "dArea" << dArea << endl;
	current_cell.store_vars(emtensor.rho, c_rho);
	current_cell.store_vars(emtensor.Si[0], c_Xmom);
        current_cell.store_vars(Stress, c_Stress);
	current_cell.store_vars(dArea, c_dArea);
    }
};

#endif /* FIXEDBGSTRESS_HPP_ */

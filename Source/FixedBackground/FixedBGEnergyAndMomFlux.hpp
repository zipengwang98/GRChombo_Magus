/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGENERGYANDMOMFLUX_HPP_
#define FIXEDBGENERGYANDMOMFLUX_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the
//! grid
template <class matter_t, class background_t> class FixedBGEnergyAndMomFlux
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    // Inherit the variable definitions from CCZ4RHS + matter_t
    template <class data_t>
    struct Vars : public CCZ4Vars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

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
    FixedBGEnergyAndMomFlux(matter_t a_matter, background_t a_background,
                            double a_dx,
                            std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        // get the metric vars from the background
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        Vars<data_t> metric_vars;
        //m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        using namespace CoordinateTransformations;
        //	const auto gamma = metric_vars.gamma;

        Tensor<2,data_t> h_UU, gamma;
        FOR(i,j){
            gamma[i][j] = vars.h[i][j] / vars.chi;
        }
        const auto gamma_UU = compute_inverse_sym(gamma);
        const auto chi = vars.chi;
        FOR(i,j){
            h_UU[i][j] = gamma_UU[i][j] / chi;
        }

        const auto chris_phys =
            compute_christoffel(d1.h, h_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, d1, h_UU, chris_phys.ULL);
        const auto lapse = metric_vars.lapse;
        const auto shift = metric_vars.shift;
        const data_t det_gamma = pow(metric_vars.chi, -3.0);
        Tensor<2, data_t> spherical_gamma = cartesian_to_spherical_LL(
            gamma, coords.x, coords.y, coords.z);
        data_t dArea = area_element_sphere(spherical_gamma);

        const data_t R = coords.get_radius();
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;

        // the unit vector in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        // Normalise
        data_t si_norm = 0.0;
        FOR2(i, j) { si_norm += gamma_UU[i][j] * si_L[i] * si_L[j]; }

        FOR1(i) { si_L[i] = si_L[i] / sqrt(si_norm); }

        Tensor<1,data_t> Mdot ;
        FOR1(i){Mdot[i] = 0;}
        FOR1(i){
            FOR1(j)
            {
                Mdot[i] += -metric_vars.shift[j] * si_L[j] * emtensor.Si[i];
                FOR1(k)
                {
                    Mdot[i] += metric_vars.lapse * gamma_UU[j][k] *
                            emtensor.Sij[i][k] * si_L[j];
                }
            }
        }
        // dArea is the integration surface element; Divide by r2sintheta,
        // as that's accounted for in the SprericalExtraction
        FOR(i){
            Mdot[i] *= dArea / r2sintheta;
        }
        

        data_t Edot = 0.0;
        FOR1(i)
        {
            Edot += lapse * si_L[i] * emtensor.rho * shift[i];
            FOR1(j)
            {
                Edot += -si_L[i] * emtensor.Si[j] *
                        (shift[i] * shift[j] + lapse * lapse * gamma_UU[i][j]);
                FOR1(k)
                {
                    Edot += si_L[i] * lapse * gamma_UU[i][j] * shift[k] *
                            emtensor.Sij[j][k];
                }
            }
        }

        Edot *= dArea / r2sintheta;

        current_cell.store_vars(Mdot[0], c_xMdot);
        current_cell.store_vars(Mdot[1], c_yMdot);
        current_cell.store_vars(Mdot[2], c_zMdot);
        current_cell.store_vars(Edot, c_Edot);
    }
};

#endif /* FIXEDBGENERGYANDMOMFLUX_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGMOMANDSOURCE_HPP_
#define FIXEDBGMOMANDSOURCE_HPP_

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
template <class matter_t, class background_t> class FixedBGMomAndSource
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
    FixedBGMomAndSource(matter_t a_matter, background_t a_background,
                        double a_dx, std::array<double, CH_SPACEDIM> a_center)
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

        //Vars<data_t> metric_vars;
        //m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        using namespace CoordinateTransformations;
        Tensor<2,data_t> h_UU, gamma;
        FOR(i,j){
            gamma[i][j] = vars.h[i][j] / vars.chi;
        }
        const auto gamma_UU = compute_inverse_sym(gamma);
        const auto chi = vars.chi;
        FOR(i,j){
            h_UU[i][j] = gamma_UU[i][j] / chi;
        }
        const auto chris =
            compute_christoffel(d1.h, h_UU);
        const auto chris_phys = compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, d1, h_UU, chris.ULL);
        const data_t det_gamma = compute_determinant_sym(gamma);

        Tensor<1,data_t> Mom;
        FOR1(i) Mom[i] = -emtensor.Si[i] * sqrt(det_gamma);

        data_t ECharge = -emtensor.rho * vars.lapse;
        FOR1(i){
            ECharge += vars.shift[i] * emtensor.Si[i];
        }
        ECharge = ECharge * sqrt(det_gamma);


        Tensor<1,data_t> Source;
        FOR1(i) Source[i] = -emtensor.rho * d1.lapse[i];
        

        FOR1(l){
            FOR1(i)
            {
                Source[l] += emtensor.Si[i] * d1.shift[i][l];
                FOR2(j, k)
                {
                    Source[l] += vars.lapse * gamma_UU[i][k] *
                            emtensor.Sij[k][j] * chris_phys[j][i][l];
                }
            }
        }

        FOR1(i) Source[i] = Source[i] * sqrt(det_gamma);


        //current_cell.store_vars(emtensor.rho, c_rho);
        current_cell.store_vars(ECharge, c_rho);
        current_cell.store_vars(Source[0], c_xSource);
        current_cell.store_vars(Source[1], c_ySource);
        current_cell.store_vars(Source[2], c_zSource);
        current_cell.store_vars(Mom[0], c_xMom);
        current_cell.store_vars(Mom[1], c_yMom);
        current_cell.store_vars(Mom[2], c_zMom);

    }
};

#endif /* FIXEDBGMOMANDSOURCE_HPP_ */

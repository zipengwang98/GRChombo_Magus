/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionDiagnostics
{
  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;
    const double m_inner_r;
    const double m_outer_r;
    const int m_loop_num;

  public:
    ExcisionDiagnostics(const double a_dx,
                        const std::array<double, CH_SPACEDIM> a_center,
                        background_t a_background, const double a_inner_r,
                        const double a_outer_r, const int a_loop_num)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_background(a_background), m_inner_r(a_inner_r), m_outer_r(a_outer_r),
          m_loop_num(a_loop_num)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        if (coords.get_radius() < m_inner_r || coords.get_radius() > m_outer_r)
        {
            if (m_loop_num == 0 ){
              current_cell.store_vars(0.0, c_xSource);
              current_cell.store_vars(0.0, c_ySource);
              current_cell.store_vars(0.0, c_xMom);
              current_cell.store_vars(0.0, c_yMom);
              current_cell.store_vars(0.0, c_rho);

            }
            else if (m_loop_num == 1 ){
              current_cell.store_vars(0.0, c_xSource_1);
              current_cell.store_vars(0.0, c_ySource_1);
              current_cell.store_vars(0.0, c_xMom_1);
              current_cell.store_vars(0.0, c_yMom_1);
              current_cell.store_vars(0.0, c_rho_1);
            }
            else if (m_loop_num == 2 ){
              current_cell.store_vars(0.0, c_xSource_2);
              current_cell.store_vars(0.0, c_ySource_2);
              current_cell.store_vars(0.0, c_xMom_2);
              current_cell.store_vars(0.0, c_yMom_2);
              current_cell.store_vars(0.0, c_rho_2);
            }
        }
    }
};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */

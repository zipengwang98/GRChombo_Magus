/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONEVOLUTION_HPP_
#define EXCISIONEVOLUTION_HPP_

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
template <class matter_t> class ExcisionEvolution
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;

  public:

    BoostedIsotropicKerr::params_t m_params;

    ExcisionEvolution(const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      BoostedIsotropicKerr::params_t a_params)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center), m_params(a_params)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);

        const double velocity = m_params.velocity;
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double boost =
            pow(1.0 - velocity * velocity, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        const double x_p = coords.x * boost;
        const double y = coords.y;
        const double z = coords.z;

        // the coordinate radius (boosted)
        const double r2 = x_p * x_p + y * y + z * z;

        // compare this to horizon in isotropic coords
        double rp = M + sqrt(M * M - a * a);
        const double r_horizon = rp / 4.0;

        double horizon_distance = sqrt(r2) / r_horizon;
        if (horizon_distance < 0.5)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;
            VarsTools::assign(vars, 0.0);
            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
        } // else do nothing
    }
};

#endif /* EXCISIONEVOLUTION_HPP_ */

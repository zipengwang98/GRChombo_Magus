/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPOTENTIAL_HPP_
#define COMPLEXPOTENTIAL_HPP_

#include "simd.hpp"

class ComplexPotential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    ComplexPotential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_re,
                           data_t &dVdphi_im, const vars_t<data_t> &vars) const
    {
        const double m = m_params.scalar_mass;
        // The potential value at phi
        // 1/2 m^2 phi^2
        V_of_phi = 0.5 * m * m * vars.phi_Re * vars.phi_Re +
	           0.5 * m * m * vars.phi_Im * vars.phi_Im;

        // The potential gradient at phi
        // m^2 phi
        dVdphi_re = m * m * vars.phi_Re;
	dVdphi_im = m * m * vars.phi_Im;
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */

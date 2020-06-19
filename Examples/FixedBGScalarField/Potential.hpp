/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
protected:
  const double m_mu;

public:
    //! The constructor
    Potential(const double a_mu) : m_mu(a_mu) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
      //        const double m = m_params.scalar_mass;
        // The potential value at phi
        // 1/2 m^2 phi^2
        V_of_phi = 0.5 * m_mu * m_mu * vars.phi * vars.phi;

        // The potential gradient at phi
        // m^2 phi
        dVdphi = m_mu * m_mu * vars.phi;
	//pout()<< "Scalar mass in potential" << m_mu <<endl;
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_phi_Re, // matter field added
    c_phi_Im,
    c_Pi_Re,  //(minus) conjugate momentum
    c_Pi_Im,
    c_chi,
    c_rho,
    c_Xmom,
    c_Stress,
    c_dArea,

    NUM_VARS
};

namespace UserVariables
{
  static constexpr char const *variable_names[NUM_VARS] = {"phi_Re", "phi_Im", "Pi_Re", "Pi_Im", "chi",
							   "rho", "Xmom", "Stress", "dArea"};
}

#endif /* USERVARIABLES_HPP */

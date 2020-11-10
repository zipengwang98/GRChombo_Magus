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
    c_Source,
    c_S1,
    c_S2,
    c_S3,
    c_d1lapse,
    c_sqdetgamma,
    c_Xmom,
    c_Stress,
    c_dArea,

    NUM_VARS
};

namespace UserVariables
{
  static constexpr char const *variable_names[NUM_VARS] = {"phi_Re", "phi_Im", "Pi_Re", "Pi_Im", "chi",
							   "rho", "Source", "S1", "S2", "S3","d1lapse" ,
							   "sqdetgamma","Xmom", "Stress", "dArea"};
}

#endif /* USERVARIABLES_HPP */

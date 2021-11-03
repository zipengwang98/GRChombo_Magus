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
    c_rho,
    c_Source,
    c_Xmom,
    c_Stress,

    NUM_VARS
}; // c_phi_Re, c_phi_Im, c_rho, c_rhofull, c_Source, c_S1, c_S2, c_S3, c_Xmom, c_Stress, c_dArea}
enum
  {
    c_empty1,
    c_empty2,
    c_empty3,
    c_empty4,
    c_empty5,
    c_empty6,
    c_Pi_Re, 
    c_Pi_Im,
    c_chi,

    NUM_VARS_CHK
  };
namespace UserVariables
{
  static constexpr char const *variable_names[NUM_VARS] = {"phi_Re", "phi_Im",
							   "rho", "Source",
							   "Xmom", "Stress"};
}

#endif /* USERVARIABLES_HPP */

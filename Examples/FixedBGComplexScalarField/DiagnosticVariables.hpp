/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_rho,
    c_xSource,
    c_ySource,
    c_zSource,
    c_xMom,
    c_yMom,
    c_zMom,
    c_xMdot, // Momentum flux
    c_yMdot, // Momentum flux
    c_zMdot, // Momentum flux
    c_Edot, // Energy flux

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "rho", "xSource", "ySource", "zSource", "xMom", "yMom", "zMom", "xMomFlux", "yMomFlux", "zMomflux", "Edot"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

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
    c_Source,
    c_xMom,
    c_xMdot, // Momentum flux
    c_yMdot, // Momentum flux
    c_zMdot, // Momentum flux
    c_Edot, // Energy flux

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "rho", "Source", "xMom", "xMomFlux", "yMomFlux", "zMomflux", "Edot"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

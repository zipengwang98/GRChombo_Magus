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
    c_xMom,
    c_yMom,
    c_rho_1,
    c_xSource_1,
    c_ySource_1,
    c_xMom_1,
    c_yMom_1,
    c_rho_2,
    c_xSource_2,
    c_ySource_2,
    c_xMom_2,
    c_yMom_2,
    c_xMdot, // Momentum flux
    c_yMdot, // Momentum flux
    c_zMdot, // Momentum flux
    c_Edot, // Energy flux

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "rho", 
    "xSource", "ySource","xMom", "yMom",
    "xSource1", "ySource1","xMom1", "yMom1",
    "xSource2", "ySource2","xMom2", "yMom2",
    "xMomFlux", "yMomFlux", "zMomflux", "Edot"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGSCALARFIELD_HPP_)
#error "This file should only be included through FixedBGScalarField.hpp"
#endif

#ifndef FIXEDBGSCALARFIELD_IMPL_HPP_
#define FIXEDBGSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // Copy the field vars into SFObject
    SFObject<data_t> vars_sf;
    vars_sf.phi = vars.phi;
    vars_sf.Pi = vars.Pi;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, metric_vars, vars_sf, d1.phi, gamma_UU,
                            chris_phys_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += -metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const SFObject<data_t> &vars_sf,
    const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars_sf.Pi * vars_sf.Pi;
    FOR2(i, j) { Vt += gamma_UU[i][j] * d1_phi[i] * d1_phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            -0.5 * metric_vars.gamma[i][j] * Vt + d1_phi[i] * d1_phi[j];
    }

    // S = Tr_S_ij
    out.S = TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i) { out.Si[i] = -d1_phi[i] * vars_sf.Pi; }

    // rho = n^a n^b T_ab
    out.rho = vars_sf.Pi * vars_sf.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGScalarField<potential_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    SFObject<data_t> rhs_sf;
    // advection terms
    SFObject<data_t> advec_sf;
    advec_sf.phi = advec.phi;
    advec_sf.Pi = advec.Pi;
    // the vars
    SFObject<data_t> vars_sf;
    vars_sf.phi = vars.phi;
    vars_sf.Pi = vars.Pi;

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(rhs_sf, vars, metric_vars, vars_sf, d1, d1.phi,
                              d2.phi, advec_sf);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term
    total_rhs.phi = rhs_sf.phi;
    total_rhs.Pi = rhs_sf.Pi - metric_vars.lapse * dVdphi;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGScalarField<potential_t>::matter_rhs_excl_potential(
    SFObject<data_t> &rhs_sf, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const SFObject<data_t> &vars_sf,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<1, data_t> &d1_phi,
    const Tensor<2, data_t> &d2_phi, const SFObject<data_t> &advec_sf)
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_sf.phi = metric_vars.lapse * vars_sf.Pi + advec_sf.phi;
    rhs_sf.Pi = metric_vars.lapse * metric_vars.K * vars_sf.Pi + advec_sf.Pi;

    FOR2(i, j)
    {
        rhs_sf.Pi += gamma_UU[i][j] * (metric_vars.lapse * d2_phi[i][j] +
                                       metric_vars.d1_lapse[i] * d1_phi[j]);
        FOR1(k)
        {
            rhs_sf.Pi += -metric_vars.lapse * gamma_UU[i][j] *
                         chris_phys.ULL[k][i][j] * d1_phi[k];
        }
    }
}

#endif /* FIXEDBGSCALARFIELD_IMPL_HPP_ */

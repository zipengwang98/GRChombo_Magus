/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGCOMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through FixedBGComplexScalarField.hpp"
#endif

#ifndef FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_
#define FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGComplexScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
  {
    // The total em tensor
    emtensor_t<data_t> out;

    // calcuate the values excluding the potential terms

    // Useful quantity Vt
    data_t Vt = -vars.Pi_Re * vars.Pi_Re - vars.Pi_Im * vars.Pi_Im;
    FOR2(i, j)
      {
        Vt += gamma_UU[i][j] * d1.phi_Re[i] * d1.phi_Re[j] +
	      gamma_UU[i][j] * d1.phi_Im[i] * d1.phi_Im[j];
      }
    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
      {
        out.Sij[i][j] =
            -0.5 * metric_vars.gamma[i][j] * Vt + d1.phi_Re[i] * d1.phi_Re[j] -
       	     0.5 * metric_vars.gamma[i][j] * Vt + d1.phi_Im[i] * d1.phi_Im[j];
      }

    // S = Tr_S_ij
    out.S = TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR1(i)
    {
        out.Si[i] = -d1.phi_Re[i] * vars.Pi_Re - d1.phi_Im[i] * vars.Pi_Im;
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi_Re * vars.Pi_Re + vars.Pi_Im * vars.Pi_Im + 0.5 * Vt;

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    // calculate total emtensor including potential terms
    out.rho = out.rho + V_of_phi;
    out.S = out.S - 3.0 * V_of_phi;
    FOR2(i, j){ out.Sij[i][j] += -metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void FixedBGComplexScalarField<potential_t>::matter_rhs(
    vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first add the terms excluding the potential
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    total_rhs.phi_Re = metric_vars.lapse * vars.Pi_Re + advec.phi_Re;
    total_rhs.Pi_Re = metric_vars.lapse * metric_vars.K * vars.Pi_Re + advec.Pi_Re;
    total_rhs.phi_Im = metric_vars.lapse * vars.Pi_Im + advec.phi_Im;
    total_rhs.Pi_Im = metric_vars.lapse * metric_vars.K * vars.Pi_Im + advec.Pi_Im;

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        total_rhs.Pi_Re +=
            gamma_UU[i][j] * (metric_vars.lapse * d2.phi_Re[i][j] +
                              metric_vars.d1_lapse[i] * d1.phi_Re[j]);
        total_rhs.Pi_Im +=
            gamma_UU[i][j] * (metric_vars.lapse * d2.phi_Im[i][j] +
                              metric_vars.d1_lapse[i] * d1.phi_Im[j]);
        FOR1(k)
        {
            total_rhs.Pi_Re += -metric_vars.lapse * gamma_UU[i][j] *
                               chris_phys.ULL[k][i][j] * d1.phi_Re[k];
            total_rhs.Pi_Im += -metric_vars.lapse * gamma_UU[i][j] *
                               chris_phys.ULL[k][i][j] * d1.phi_Im[k];
        }
    }

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    // adjust RHS for Pi for this
    total_rhs.Pi_Re += -metric_vars.lapse * dVdphi_Re;
    total_rhs.Pi_Im += -metric_vars.lapse * dVdphi_Im;
}

#endif /* FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_ */

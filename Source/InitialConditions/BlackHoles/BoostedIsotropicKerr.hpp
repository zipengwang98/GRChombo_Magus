/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDISOTROPICKERR_HPP_
#define BOOSTEDISOTROPICKERR_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Class which computes the initial conditions for a boosted BH
// in isotropic Schwazschild coords
// Undoosted elem. ds^2 = -Adt + B dx^2
// Lorentz boost x = \gamma (x' - vt'), t = \gamma (t' - v*x')
// Addisional shift to fix coordinates to the BH
// x' = \tilde x + v \tilde t, t' = \tilde t
// BH appears fixed, but we're in the rest frame of the SF!

class BoostedIsotropicKerr
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        double spin = 0.0;                      //!<< The spin of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double velocity = 0.0; //!< The boost velocity in the x direction
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    BoostedIsotropicKerr(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
        // check this boost param is sensible
        if ((m_params.velocity > 1.0) || (m_params.velocity < -1.0))
        {
            MayDay::Error("The boost velocity parameter must be in the range "
                          "-1.0 < velocity < 1.0");
        }
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error("The spin parameter must be in the range "
                          "-M < a< M");
        }
    }

    template <class data_t>
    data_t get_func(double M, double a, double V, data_t x, double y, double z, int index) const;

    template <class data_t>
    data_t get_d_func(double M, double a, double V, data_t x, double y, double z, int index, int d_index) const;

    template <class data_t>
    data_t get_alpha(double M, double a, double V, data_t x, double y, double z) const; 

    template<class data_t>
    Tensor<1,data_t> get_beta(double M, double a, double V, data_t x, double y, double z) const;

    template<class data_t>
    Tensor<2,data_t> get_gamma(double M, double a, double V, data_t x, double y, double z) const;

    template<class data_t>
    Tensor<1,data_t> get_d1alpha(double M, double a, double V, data_t x, double y, double z) const;

    template<class data_t>
    Tensor<2,data_t> get_d1beta(double M, double a, double V, data_t x, double y, double z) const;

    template<class data_t>
    Tensor<3,data_t> get_d1gamma(double M, double a, double V, data_t x, double y, double z) const;

    template<class data_t>
    const Tensor<2,data_t> get_K(double M, double a, double V, data_t x, double y, double z) const;

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    /// Schwarzschild boosted solution as above
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const;


  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const;
    
};
#include "BoostedIsotropicKerr.impl.hpp"
#endif /* BOOSTEDISOTROPICKERR_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARGAUSSIAN_HPP_
#define SCALARGAUSSIAN_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a scalar field given params for initial
//! matter config
class ScalarGaussian
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
    };

    //! The constructor
    ScalarGaussian(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        ScalarField<>::Vars<data_t> vars;
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t r = coords.get_radius();

        // set the field vars
        vars.phi =
            m_params.amplitude * exp(-r * r / m_params.width / m_params.width);
        vars.Pi = 0;

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* SCALARGAUSSIAN_HPP_ */

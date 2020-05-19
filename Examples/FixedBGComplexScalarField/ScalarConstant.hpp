/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARCONSTANT_HPP_
#define SCALARCONSTANT_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
class ScalarConstant
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
      //std::array<double, CH_SPACEDIM>
	//        center;   //!< Centre of perturbation in initial SF bubble
	//    double width; //!< Width of bump in initial SF bubble
    };

    //! The constructor
    ScalarConstant(params_t a_params)
        : m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        FixedBGComplexScalarField<>::Vars<data_t> vars;
        //Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        //data_t r = coords.get_radius();
	const double omega = 1.;
        // set the field vars
        vars.phi_Re = m_params.amplitude; //* exp(-r * r / m_params.width / m_params.width);
	vars.phi_Im = 0;
        vars.Pi_Re = 0;
	vars.Pi_Im = m_params.amplitude * omega;

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }

  protected:
  //double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* SCALARCONSTANT_HPP_ */

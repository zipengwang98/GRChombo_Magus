/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CUSTOMEXTRACTION_HPP_
#define CUSTOMEXTRACTION_HPP_
#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParameters.hpp"
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp"
#include <fstream>
#include <iostream>

//!  The class allows extraction of the values of components at
//!  specified points
class CustomExtraction
{
private:
  //! Params for extraction
  const int m_comp;
  const int m_num_points;
  const double m_L;
  const std::array<double, CH_SPACEDIM> m_center;
  const double m_dt;
  const double m_time;

public:
  //! The constructor
  CustomExtraction(int a_comp, int a_num_points, double a_L,
		   std::array<double, CH_SPACEDIM> a_center, double a_dt,
		   double a_time)
    : m_comp(a_comp), m_num_points(a_num_points), m_center(a_center),
      m_L(a_L), m_dt(a_dt), m_time(a_time)
  {
  }

  //! Destructor
  ~CustomExtraction() {}

  //! Execute the query
  void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator,
		     std::string a_filename) const;

private:
  //! Write out the result of the extraction in phi and theta at each timestep
  //! for each extraction radius
  void write_extraction(std::string a_filename,
			const std::vector<double> &interp_x,
			const std::vector<double> &interp_y,
			const std::vector<double> &interp_z,
			const std::vector<double> &var_data) const;
};

#include "CustomExtraction.impl.hpp"

#endif /* CUSTOMEXTRACTION_HPP_ */ 

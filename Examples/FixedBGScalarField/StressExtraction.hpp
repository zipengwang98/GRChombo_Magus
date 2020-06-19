/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STRESSEXTRACTION_HPP_
#define STRESSEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParameters.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp" // Needs c_Stress and c_dArea
//!  The class allows extraction of the values of the Stress scalar components on
//!  spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the Stress
   over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class StressExtraction
{
  private:
    //! Params for extraction
    const extraction_params_t m_params;
    const int m_Stress = c_Stress;
    const int m_dArea = c_dArea;
  //const int m_rho = c_rho;
    const double m_dt;
    const double m_time;
    const bool m_first_step;
    const double m_restart_time;
    const int m_num_points; // number of points per extraction radius
    const double m_dphi;
    const double m_dtheta;

  public:
    //! The constructor
    StressExtraction(extraction_params_t a_params, double a_dt, double a_time,
                   bool a_first_step, double a_restart_time = 0.0)
        : m_params(a_params), m_dt(a_dt), m_time(a_time),
          m_first_step(a_first_step), m_restart_time(a_restart_time),
          m_num_points(m_params.num_points_phi * m_params.num_points_theta),
          m_dphi(M_PI / m_params.num_points_phi),
          m_dtheta(M_PI /2.0 / m_params.num_points_theta)
    {
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    StressExtraction(extraction_params_t a_params, double a_dt, double a_time,
                   double a_restart_time = 0.0)
        : StressExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }


    //! Destructor
    ~StressExtraction() {}

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator) const;

  private:
    //! integrate over a spherical shell with given harmonics for each
    //! extraction radius and normalise by multiplying by radius
    std::pair<std::vector<double>, std::vector<double>>
    integrate_surface(int es, int el, int em,
                      const std::vector<double> a_Stress,
                      const std::vector<double> a_dArea) const;

    //! Write out calculated values of integral for each extraction radius
    void write_integral(const std::vector<double> a_integral_Stress_1,
                        const std::vector<double> a_integral_Stress_2,
                        std::string a_filename) const;

    //! Write out the result of the extraction in phi and theta at each timestep
    //! for each extraction radius
    void write_extraction(std::string a_file_prefix,
                          const std::vector<double> a_Stress,
                          const std::vector<double> a_dArea) const;
};

#include "StressExtraction.impl.hpp"

#endif /* STRESSEXTRACTION_HPP_ */

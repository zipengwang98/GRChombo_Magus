/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(STRESSEXTRACTION_HPP_)
#error "This file should only be included through StressExtraction.hpp"
#endif

#ifndef STRESSEXTRACTION_IMPL_HPP_
#define STRESSEXTRACTION_IMPL_HPP_

//! Set up and execute the interpolation query
inline void StressExtraction::execute_query(
    AMRInterpolator<Lagrange<4>> *a_interpolator) const
{
    CH_TIME("StressExtraction::execute_query");
    if (a_interpolator == nullptr)
    {
        MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }

    std::vector<double> interp_Stress(m_num_points *
                                       m_params.num_extraction_radii);
    std::vector<double> interp_Momx(m_num_points *
                                       m_params.num_extraction_radii);
    std::vector<double> interp_x(m_num_points * m_params.num_extraction_radii);
    std::vector<double> interp_y(m_num_points * m_params.num_extraction_radii);
    std::vector<double> interp_z(m_num_points * m_params.num_extraction_radii);

    // Work out the coordinates
    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = (iphi + 0.5) * m_dphi;
            interp_x[iradius * m_num_points + idx] =
                m_params.extraction_center[0] +
                m_params.extraction_radii[iradius] * sin(theta) * cos(phi);
            interp_y[iradius * m_num_points + idx] =
                m_params.extraction_center[1] +
                m_params.extraction_radii[iradius] * sin(theta) * sin(phi);
            interp_z[iradius * m_num_points + idx] =
                m_params.extraction_center[2] +
                m_params.extraction_radii[iradius] * cos(theta);
	}
    }
    // set up the query
    InterpolationQuery query(m_num_points * m_params.num_extraction_radii);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(m_Stress, interp_Stress.data())
        .addComp(m_Momx, interp_Momx.data());

    // submit the query
    a_interpolator->interp(query);

    for (int imode = 0; imode < m_params.num_modes; ++imode)
    {
        const std::pair<int, int> &mode = m_params.modes[imode];
        auto integral = integrate_surface(-2, mode.first, mode.second,
                                          interp_Stress, interp_Momx);
        std::string integral_filename = "Stress_integral_" +
                                        std::to_string(mode.first) +
                                        std::to_string(mode.second);
        write_integral(integral.first, integral.second, integral_filename);
    }

    if (m_params.write_extraction)
    {
        write_extraction("ExtractionOut_", interp_Stress, interp_Momx);
    }
}

//! integrate over a spherical shell with given harmonics for each extraction
//! radius and normalise by multiplying by radius
inline std::pair<std::vector<double>, std::vector<double>>
StressExtraction::integrate_surface(int es, int el, int em,
                                  const std::vector<double> a_Stress,
				  const std::vector<double> a_Momx) const
{
    CH_TIME("StressExtraction::integrate_surface");
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    std::vector<double> integral_Stress(m_params.num_extraction_radii, 0.);
    //std::vector<double> inner_integral_Momx(m_params.num_extraction_radii, 0.);
    std::vector<double> integral_Momx(m_params.num_extraction_radii, 0.);
    double inner_integral_Momx = 0.;  
    // only rank 0 does the integral, but use OMP threads if available
    if (rank == 0)
    {
        // integrate the values over the sphere (normalised by r^2)
        // for each radius
        // assumes spacings constant, uses trapezium rule for phi and rectangles
        // for theta  note we don't have to fudge the end points for phi because
        // the function is periodic  and so the last point (implied but not part
        // of vector) is equal to the first point
#ifdef _OPENMP
#if __GNUC__ > 8
#define OPENMP_CONST_SHARED shared(a_Stress, a_Momx)
#else
#define OPENMP_CONST_SHARED
#endif
#pragma omp parallel for collapse(2) default(none)                             \
  shared(es, el, em, integral_Stress, integral_Momx) OPENMP_CONST_SHARED
#undef OPENMP_CONST_SHARED
#endif
      for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
	{
	  //	  double inner_integral_Momx = 0.;
	  for (int iphi = 0; iphi < m_params.num_points_phi; ++iphi)
            {
	      double inner_integral_Momx = 0.;
	      double phi = (iphi + 0.5) * m_dphi;
	      double inner_integral_Stress = 0.;
	      double f_theta_phi_r_Momx = 0.;
	      for (int itheta = 0; itheta < m_params.num_points_theta;
		   itheta++)
                {
		  using namespace SphericalHarmonics;
		  double theta = (itheta + 0.5) * m_dtheta;
		  int idx = iradius * m_num_points +
		    itheta * m_params.num_points_phi + iphi;
		  double x = m_params.extraction_radii[iradius] * sin(theta) *
		    cos(phi);
		  double y = m_params.extraction_radii[iradius] * sin(theta) *
		    sin(phi);
		  double z = m_params.extraction_radii[iradius] * cos(theta);
		  
		  double integrand_Stress = a_Stress[idx];
		  double integrand_Momx = a_Momx[idx];
		  double r2sintheta =  m_params.extraction_radii[iradius] *  m_params.extraction_radii[iradius] * sin(theta);
		  double f_theta_phi_Stress = integrand_Stress * r2sintheta;
		  double f_theta_phi_Momx = integrand_Momx * r2sintheta;
		  inner_integral_Stress += m_dtheta * f_theta_phi_Stress;
		  f_theta_phi_r_Momx += m_dtheta * f_theta_phi_Momx;
                }
#ifdef _OPENMP
#pragma omp atomic
#endif
	      integral_Stress[iradius] += m_dphi * inner_integral_Stress;
	      inner_integral_Momx += m_dphi * f_theta_phi_r_Momx;
#ifdef _OPENMP
#pragma omp atomic
#endif
	      integral_Momx[iradius] += m_dx * inner_integral_Momx;
	    }
        }
    }
    return std::make_pair(integral_Stress, integral_Momx);
}

//! Write out calculated value of integral for each extraction radius
inline void
StressExtraction::write_integral(const std::vector<double> a_integral_Stress,
                               const std::vector<double> a_integral_Momx,
                               std::string a_filename) const
{
    CH_TIME("StressExtraction::write_integral");
    // open file for writing
    SmallDataIO integral_file(a_filename, m_dt, m_time, m_restart_time,
                              SmallDataIO::APPEND, m_first_step);

    // remove any duplicate data if this is a restart
    integral_file.remove_duplicate_time_data();

    // need to write headers if this is the first timestep
    if (m_first_step)
    {
        // make header strings
        std::vector<std::string> header1_strings(2 *
                                                 m_params.num_extraction_radii);
        std::vector<std::string> header2_strings(2 *
                                                 m_params.num_extraction_radii);
        for (int iintegral = 0; iintegral < 2 * m_params.num_extraction_radii;
             iintegral += 2)
        {
            int iintegral1 = iintegral + 1;
            int iradius = iintegral / 2;
            header1_strings[iintegral] = "integral Stress";
            header1_strings[iintegral1] = "integral Mom_x";
            header2_strings[iintegral1] = header2_strings[iintegral] =
                std::to_string(m_params.extraction_radii[iradius]);
        }

        // write headers
        integral_file.write_header_line(header1_strings);
        integral_file.write_header_line(header2_strings, "r = ");
    }

    // make vector of data for writing
    std::vector<double> data_for_writing(2 * m_params.num_extraction_radii);
    for (int iintegral = 0; iintegral < 2 * m_params.num_extraction_radii;
         iintegral += 2)
    {
        int iintegral1 = iintegral + 1;
        int iradius = iintegral / 2;
        data_for_writing[iintegral] = a_integral_Stress[iradius];
        data_for_writing[iintegral1] = a_integral_Momx[iradius];
    }

    // write data
    integral_file.write_time_data_line(data_for_writing);
}

//! Write out the result of the extraction in phi and theta at each timestep for
//! each extraction radius
inline void
StressExtraction::write_extraction(std::string a_file_prefix,
                                   const std::vector<double> a_Stress,
				   const std::vector<double> a_Momx) const
{
    CH_TIME("StressExtraction::write_extraction");
    SmallDataIO extraction_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW, m_first_step);

    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        // Write headers
        std::vector<std::string> header1_strings = {
            "time = " + std::to_string(m_time) + ",",
            "r = " + std::to_string(m_params.extraction_radii[iradius])};
        extraction_file.write_header_line(header1_strings, "");
        std::vector<std::string> components = {
            UserVariables::variable_names[m_Stress],
            UserVariables::variable_names[m_Momx]};
        std::vector<std::string> coords = {"theta", "phi"};
        extraction_file.write_header_line(components, coords);

        // Now the data
        for (int idx = iradius * m_num_points;
             idx < (iradius + 1) * m_num_points; ++idx)
        {
            int itheta =
                (idx - iradius * m_num_points) / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = (iphi + 0.5) * m_dphi;

            extraction_file.write_data_line({a_Stress[idx], a_Momx[idx]},
                                            {theta, phi});
        }
        extraction_file.line_break();
    }
}

#endif /* STRESSEXTRACTION_IMPL_HPP_ */

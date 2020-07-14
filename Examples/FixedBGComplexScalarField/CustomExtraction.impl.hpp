/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CUSTOMEXTRACTION_HPP_)
#error "This file should only be included through CustomExtraction.hpp"
#endif

#ifndef CUSTOMEXTRACTION_IMPL_HPP_
#define CUSTOMEXTRACTION_IMPL_HPP_

//! Set up and execute the interpolation query
inline void CustomExtraction::execute_query(
					    AMRInterpolator<Lagrange<4>> *a_interpolator, std::string a_filename) const
{
  CH_TIME("CustomExtraction::execute_query");
  if (a_interpolator == nullptr)
    {
      MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }
  std::vector<double> interp_var_data(m_num_points);
  std::vector<double> interp_x(m_num_points);
  std::vector<double> interp_y(m_num_points);
  std::vector<double> interp_z(m_num_points);

  // Work out the coordinates
  // go out radially, focussed near centre, out to L/2
  for (int idx = 0; idx < m_num_points; ++idx)
    {
        interp_x[idx] =
            m_center[0] + 1.0 +
	  pow(double(idx) / double(m_num_points), 3.0) * 0.5 * m_L;
        interp_y[idx] =
            m_center[1] + 1.0 +
	  pow(double(idx) / double(m_num_points), 3.0) * 0.5 * m_L;
        interp_z[idx] =
            m_center[2] + 1.0 +
	  pow(double(idx) / double(m_num_points), 3.0) * 0.5 * m_L;
    }

  // set up the query
  InterpolationQuery query(m_num_points);
  query.setCoords(0, interp_x.data())
    .setCoords(1, interp_y.data())
    .setCoords(2, interp_z.data())
    .addComp(m_comp, interp_var_data.data());

  // submit the query
  a_interpolator->interp(query);

  // This generates the output
  write_extraction(a_filename, interp_x, interp_y, interp_z, interp_var_data);
}

//! Write out the result of the extraction in phi and theta at each timestep for
//! each extraction radius
inline void CustomExtraction::write_extraction(
					       std::string a_filename, const std::vector<double> &interp_x,
					       const std::vector<double> &interp_y, const std::vector<double> &interp_z,
					       const std::vector<double> &var_data) const
{
  CH_TIME("CustomExtraction::write_extraction");
  int rank;
#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
  rank = 0;
#endif
  // only rank 0 does the write out
  if (rank == 0)
    {
      // set up component names
      std::string comp_str = UserVariables::variable_names[m_comp];

      // write out extraction data to same file at each step
      std::ofstream outfile;
      if (m_time == m_dt)
        {
	  outfile.open(a_filename);
        }
      else
        {
	  outfile.open(a_filename, std::ios_base::app);
        }
      if (!outfile)
        {
	  MayDay::Error("CustomExtraction::write_extraction: error opening "
			"output file");
        }

      // header data if first timestep
      if (m_time == m_dt)
        {
	  outfile << "# Component = " << comp_str << "\n";
	  outfile << "#" << std::setw(19) << "time, r = ";
	  for (int idx = 0; idx < m_num_points; idx++)
            {
	      double r = sqrt(interp_x[idx] * interp_x[idx] +
                                interp_y[idx] * interp_y[idx] +
			      interp_z[idx] * interp_z[idx]);
	      outfile << std::setw(20) << r;
            }
	  outfile << "\n";
        }

      outfile << std::setw(20) << m_time;
      for (int idx = 0; idx < m_num_points; idx++)
        {
	  outfile << std::setw(20) << var_data[idx];
        }
      outfile << "\n";
      outfile.close();
    }
}

#endif /* CUSTOMEXTRACTION_IMPL_HPP_ */

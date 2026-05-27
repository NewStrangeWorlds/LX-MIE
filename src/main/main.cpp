/*
* This file is part of the LX-MIE code (https://github.com/exoclime).
* Copyright (C) 2017 Daniel Kitzmann
*
* LX-MIE is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* LX-MIE is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* LX-MIE directory under <license.md>. If not, see
* <http://www.gnu.org/licenses/>.
*/



#include "../mie/mie.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <cmath>
#include <stdexcept>


#include "main_auxiliary.h"



int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: lx_mie <parameter_file>\n";
    return 1;
  }

  std::string parameter_file = argv[1];

  std::string refractive_index_file;
  std::string output_file;
  double particle_radius;
  bool use_cross_sections = true;
  bool compute_phase_function = false;
  bool use_legendre_series = false;

  std::vector<double> phase_function_angles;
  std::string phase_function_output_file = "";
  size_t nb_legendre_terms = 0;



  if (readParameterFile(parameter_file, refractive_index_file, output_file, particle_radius,
                        use_cross_sections, compute_phase_function, use_legendre_series,
                        phase_function_angles, nb_legendre_terms, phase_function_output_file) == false) return 1;



  std::vector< std::complex<double> > refractive_index;
  std::vector<double> wavelengths;

  if (readRefractiveIndex(refractive_index_file, wavelengths, refractive_index) == false) return 1;



  size_t nb_wavelengths = wavelengths.size();


  std::vector<double> q_abs(nb_wavelengths, 0);
  std::vector<double> q_ext(nb_wavelengths, 0);
  std::vector<double> q_scat(nb_wavelengths, 0);
  std::vector<double> asymmetry_parameter(nb_wavelengths, 0);


  std::vector< std::vector<double> > phase_function;
  phase_function.resize(nb_wavelengths);


  try
  {
    if (!compute_phase_function)
    {
      for (size_t i = 0; i < nb_wavelengths; ++i)
      {
        double size_parameter = 2.0 * CONST_PI * particle_radius / wavelengths[i];

        auto r = lxmie::Mie(refractive_index[i], size_parameter);

        q_ext[i]              = r.q_ext;
        q_scat[i]             = r.q_sca;
        q_abs[i]              = r.q_abs;
        asymmetry_parameter[i] = r.asymmetry_parameter;
      }
    }

    if (compute_phase_function)
    {
      if (use_legendre_series)
      {
        for (size_t i = 0; i < nb_wavelengths; ++i)
        {
          double size_parameter = 2.0 * CONST_PI * particle_radius / wavelengths[i];

          phase_function[i].assign(nb_legendre_terms, 0);

          auto r = lxmie::Mie(refractive_index[i], size_parameter,
                               nb_legendre_terms, phase_function[i]);

          q_ext[i]              = r.q_ext;
          q_scat[i]             = r.q_sca;
          q_abs[i]              = r.q_abs;
          asymmetry_parameter[i] = r.asymmetry_parameter;
        }

        writeLegendreSeriesOutput(phase_function_output_file, wavelengths, phase_function);
      }
      else
      {
        for (size_t i = 0; i < nb_wavelengths; ++i)
        {
          double size_parameter = 2.0 * CONST_PI * particle_radius / wavelengths[i];

          phase_function[i].assign(phase_function_angles.size(), 0);

          std::vector< std::complex<double> > s1(phase_function_angles.size(), std::complex<double>(0, 0));
          std::vector< std::complex<double> > s2(phase_function_angles.size(), std::complex<double>(0, 0));

          auto r = lxmie::Mie(refractive_index[i], size_parameter, phase_function_angles, s1, s2);

          q_ext[i]              = r.q_ext;
          q_scat[i]             = r.q_sca;
          q_abs[i]              = r.q_abs;
          asymmetry_parameter[i] = r.asymmetry_parameter;

          phase_function[i] = lxmie::phaseFunction(s1, s2, size_parameter, q_scat[i]);
        }

        writePhaseFunctionOutput(phase_function_output_file, wavelengths, phase_function_angles, phase_function);
      }
    }
  }
  catch (const std::runtime_error& e)
  {
    std::cerr << "Mie calculation error: " << e.what() << "\n";
    return 1;
  }


  writeOutputFile(output_file, use_cross_sections, particle_radius, wavelengths, q_ext, q_abs, q_scat, asymmetry_parameter);



  std::cout << "Model finished!\n";

  return 0;
}

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



#include "main_auxiliary.h"


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <cmath>


//write results
void writeOutputFile(const std::string file_name, const bool use_cross_sections, const double particle_radius,
                     const std::vector<double> wavelengths,
                     const std::vector<double>& q_ext, const std::vector<double>& q_abs,
                     const std::vector<double>& q_scat, const std::vector<double>& asymmetry_parameter)
{
  std::fstream file(file_name.c_str(), std::ios::out);


  if (file.fail())
  {
    std::cout << "Could not create output file: " << file_name << ". You wasted a lot of energy! :.-( \n";

    return;
  }


  size_t nb_wavelengths = wavelengths.size();


  if (use_cross_sections)
  {
    double radius = particle_radius *1e-4; //convert from micron to cm to calculate cross section in cm2

    file << std::setw(16) << std::left << "#wavelengths (mu)" << "\t"
         << std::setw(16) << std::left << "size parameter" << "\t"
         << std::setw(33) << std::left << "extinction cross section (cm^2)" << "\t"
         << std::setw(33) << std::left << "scattering cross section (cm^2)" << "\t"
         << std::setw(33) << std::left << "absorption cross section (cm^2)" << "\t"
         << std::setw(33) << std::left << "single scattering albedo" << "\t"
         << std::setw(33) << std::left << "asymmetry parameter" << "\n";

    for (size_t i=0; i<nb_wavelengths; i++)
      file << std::setw(16) << std::left << std::setprecision(10) << std::scientific << wavelengths[i] << "\t"
           << std::setw(16) << std::left << 2.0 * CONST_PI * particle_radius / wavelengths[i] << "\t"
           << std::setw(33) << std::left << q_ext[i]*CONST_PI*radius*radius<< "\t"
           << std::setw(33) << std::left << q_scat[i]*CONST_PI*radius*radius << "\t"
           << std::setw(33) << std::left << q_abs[i]*CONST_PI*radius*radius << "\t"
           << std::setw(33) << std::left << q_scat[i]/q_ext[i] << "\t"
           << std::setw(33) << std::left << asymmetry_parameter[i] << "\n";
  }
  else
  {
    file << std::setw(16) << std::left << "#wavelengths (mu)" << "\t"
         << std::setw(16) << std::left << "size parameter" << "\t"
         << std::setw(16) << std::left << "q_abs" << "\t"
         << std::setw(16) << std::left << "q_scat" << "\t"
         << std::setw(16) << std::left << "q_ext" << "\t"
         << std::setw(16) << std::left << "asymmetry parameter" << "\n";

    for (size_t i=0; i<nb_wavelengths; i++)
      file << std::setprecision(10) << std::scientific << wavelengths[i] << "\t"
                                                       << 2.0 * CONST_PI * particle_radius / wavelengths[i] << "\t"
                                                       << q_abs[i] << "\t"
                                                       << q_scat[i] << "\t"
                                                       << q_ext[i] << "\t"
                                                       << asymmetry_parameter[i] << "\n";
  }


  file.close();
}


void writePhaseFunctionOutput(const std::string file_name, const std::vector<double>& wavelengths,
                              const std::vector<double>& angles, const std::vector< std::vector<double> >& phase_function)
{
  std::fstream file(file_name.c_str(), std::ios::out);


  if (file.fail())
  {
    std::cout << "Could not create phase function output file: " << file_name << ". You wasted a lot of energy! :.-( \n";

    return;
  }

  size_t nb_wavelengths = wavelengths.size();
  size_t nb_angles = angles.size();


  file << std::setw(16) << std::left << "#cos(angle)" << "\t"
       << std::setw(16) << std::left << "phase function" << "\n";

  for (size_t i=0; i<nb_angles; ++i)
  {
    file << std::setprecision(10) << std::scientific << angles[i];

    for (size_t j=0; j<nb_wavelengths; ++j)
      file << "\t" << std::setprecision(10) << std::scientific << phase_function[j][i];

    file << "\n";
  }


}





void writeLegendreSeriesOutput(const std::string file_name, const std::vector<double>& wavelengths,
                               const std::vector< std::vector<double> >& legendre_series)
{
  std::fstream file(file_name.c_str(), std::ios::out);


  if (file.fail())
  {
    std::cout << "Could not create phase function output file: " << file_name << ". You wasted a lot of energy! :.-( \n";

    return;
  }

  size_t nb_wavelengths = wavelengths.size();
  size_t nb_terms = legendre_series[0].size();


  file << std::setw(16) << std::left << "#Legendre series terms" << "\n";

  for (size_t i=0; i<nb_terms; ++i)
  {
    file << std::setprecision(6) << std::scientific << i;

    for (size_t j=0; j<nb_wavelengths; ++j)
      file << "\t" << std::setprecision(10) << std::scientific << legendre_series[j][i];

    file << "\n";
  }


}




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



#ifndef MAIN_AUXILIARY_H_
#define MAIN_AUXILIARY_H_

#include <vector>
#include <complex>
#include <fstream>



double const CONST_PI = 3.14159265358979323846;


bool readParameterFile(const std::string parameter_file, std::string &refractive_index_file, std::string &output_file, double &radius, bool& use_cross_sections,
                       bool& compute_phase_function, bool& use_legendre_series, std::vector<double>& phase_function_angles, size_t& nb_legendre_terms,
                       std::string& phase_function_output_file);

bool readPhaseFunctionParameter(std::fstream& file, const bool use_legendre_series, std::vector<double>& phase_angles, size_t& nb_legendre_terms,
                                std::string& output_file);

bool readRefractiveIndex(const std::string file_name, std::vector<double> &wavelengths, std::vector< std::complex<double> > &refractive_index);


void writeOutputFile(const std::string file_name, const bool use_cross_sections, const double particle_radius,
                     const std::vector<double> wavelengths,
                     const std::vector<double>& q_ext, const std::vector<double>& q_abs,
                     const std::vector<double>& q_scat, const std::vector<double>& asymmetry_parameter);

void writePhaseFunctionOutput(const std::string file_name, const std::vector<double>& wavelengths,
                              const std::vector<double>& angles, const std::vector< std::vector<double> >& phase_function);

void writeLegendreSeriesOutput(const std::string file_name, const std::vector<double>& wavelengths,
                               const std::vector< std::vector<double> >& legendre_series);


#endif


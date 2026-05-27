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



bool readPhaseFunctionParameter(std::fstream& file, const bool use_legendre_series, std::vector<double>& phase_angles, size_t& nb_legendre_terms,
                                std::string& output_file)
{
  std::string line;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  file >> output_file;
  if (output_file == "") return false;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  if (!use_legendre_series)
  {
    phase_angles.resize(0);


    while(std::getline(file, line))
    {
      if (file.eof()) break;

      std::stringstream  line_stream(line);

      double angle;

      line_stream >> angle;

      phase_angles.push_back(angle);
    }

    if (phase_angles.size() == 0)
      return false;
  }
  else
  {
    file >> nb_legendre_terms;

    if (nb_legendre_terms == 0)
      return false;
  }



  return true;
}



bool readParameterFile(const std::string parameter_file, std::string &refractive_index_file, std::string &output_file, double &radius, bool& use_cross_sections,
                       bool& compute_phase_function, bool& use_legendre_series, std::vector<double>& phase_function_angles, size_t& nb_legendre_terms,
                       std::string& phase_function_output_file)
{
  std::fstream file(parameter_file.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "parameter file " << parameter_file << " not found! :-( \n";
    return false;
  }

  std::string line;

  std::getline(file, line);

  file >> refractive_index_file;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  file >> radius;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  file >> output_file;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  file >> line;
  if (line == "CS")
    use_cross_sections = true;
  else
    use_cross_sections = false;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  if (!file.eof())
  {
    file >> line;

    if (line == "PF" || line == "LS")
    {
      compute_phase_function = true;

      if (line == "LS")
        use_legendre_series = true;
      else
        use_legendre_series = false;


      compute_phase_function = readPhaseFunctionParameter(file, use_legendre_series, phase_function_angles, nb_legendre_terms, phase_function_output_file);

      if (!compute_phase_function)
        std::cout << "Error reading the phase function/Legendre series parameters!\n";

    }
    else
      compute_phase_function = false;
  }


  file.close();


  std::cout << "read from parameter file: \n";
  std::cout << "refractive index: " << refractive_index_file << "\n";
  std::cout << "particle radius: " << radius << "\n";
  std::cout << "output file: " << output_file << "\n";
  std::cout << "use cross section " << use_cross_sections << "\n";
  std::cout << "compute phase function " << compute_phase_function << "\n";

  if (compute_phase_function)
  {
    std::cout << "phase function output " << phase_function_output_file << "\n";

    if (!use_legendre_series)
    {
      std::cout << "phase function angles: \n";

      for (size_t i=0; i<phase_function_angles.size(); ++i)
        std::cout << phase_function_angles[i] << "\n";
    }
    else
    {
      std::cout << "Number of terms in Legendre series: " << nb_legendre_terms << "\n";
    }

  }



  return true;
}




bool readRefractiveIndex(const std::string file_name, std::vector<double> &wavelengths, std::vector< std::complex<double> > &refractive_index)
{
  std::fstream file(file_name.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "refractive index file " << file_name << " not found! :-( \n";
    return false;
  }


  std::string line;

  //file comments
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  double wavelength, n, k;

  while (file >> wavelength >> n >> k)
  {
    if (k>0) k *= -1.;

    wavelengths.push_back(wavelength);
    refractive_index.push_back(std::complex<double>(n, k));
  }


  file.close();


  if (wavelengths.size() == 0)
    return false;


  return true;
}


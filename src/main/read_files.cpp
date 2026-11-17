/*
* This file is part of the LX-MIE code (https://github.com/newstrangeworlds).
* Copyright (C) 2026 Daniel Kitzmann
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



bool readPhaseFunctionParameter(
  std::fstream& file, 
  const bool use_legendre_series, 
  std::vector<double>& phase_angles, 
  size_t& nb_legendre_terms,
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



//Reads the particle size from the corresponding line of the parameter file
//This is either a single particle radius or one of the following size distributions:
//  lognormal <median radius> <geometric standard deviation> [number of radii]
//  gamma     <effective radius> <effective variance> [number of radii]
//  powerlaw  <minimum radius> <maximum radius> <exponent> [number of radii]
//  file      <file name>
bool readParticleSize(
  const std::string& input, 
  double& radius, 
  bool& use_size_distribution,
  lxmie::SizeDistribution& size_distribution)
{
  std::stringstream input_stream(input);

  //a single particle radius
  if (input_stream >> radius)
  {
    use_size_distribution = false;

    if (radius > 0) return true;

    std::cout << "The particle radius has to be positive!\n";
    return false;
  }

  use_size_distribution = true;

  input_stream = std::stringstream(input);

  std::string type;
  input_stream >> type;

  try
  {
    if (type == "file")
    {
      std::string file_name;

      if (input_stream >> file_name)
        return readSizeDistributionFile(file_name, size_distribution);
    }
    else if (type == "lognormal" || type == "gamma")
    {
      double parameter_1 = 0, parameter_2 = 0;
      size_t nb_points = lxmie::default_nb_distribution_points;

      if (input_stream >> parameter_1 >> parameter_2)
      {
        input_stream >> nb_points;

        if (type == "lognormal")
          size_distribution = lxmie::logNormalDistribution(parameter_1, parameter_2, nb_points);
        else
          size_distribution = lxmie::gammaDistribution(parameter_1, parameter_2, nb_points);

        return true;
      }
    }
    else if (type == "powerlaw")
    {
      double r_min = 0, r_max = 0, exponent = 0;
      size_t nb_points = lxmie::default_nb_distribution_points;

      if (input_stream >> r_min >> r_max >> exponent)
      {
        input_stream >> nb_points;

        size_distribution = lxmie::powerLawDistribution(r_min, r_max, exponent, nb_points);

        return true;
      }
    }
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << "\n";
    return false;
  }

  std::cout << "Error reading the particle size from the parameter file: " << input << "\n";

  return false;
}



//Reads a tabulated size distribution
//The file contains two columns: particle radius in micron and the distribution dN/dr
bool readSizeDistributionFile(
  const std::string file_name, 
  lxmie::SizeDistribution& size_distribution)
{
  std::fstream file(file_name.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "size distribution file " << file_name << " not found! :-( \n";
    return false;
  }


  std::vector<double> radius;
  std::vector<double> dn_dr;

  std::string line;

  while (std::getline(file, line))
  {
    const size_t first = line.find_first_not_of(" \t\r");

    if (first == std::string::npos || line[first] == '#') continue;

    std::stringstream line_stream(line);
    double r, n;

    if (!(line_stream >> r >> n))
    {
      std::cout << "Error reading the size distribution file " << file_name << ": " << line << "\n";
      return false;
    }

    radius.push_back(r);
    dn_dr.push_back(n);
  }


  file.close();


  try
  {
    size_distribution = lxmie::tabulatedDistribution(radius, dn_dr);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << "\n";
    return false;
  }


  return true;
}



bool readParameterFile(
  const std::string parameter_file, 
  std::string &refractive_index_file, 
  std::string &output_file, 
  double &radius, 
  bool& use_size_distribution,
  lxmie::SizeDistribution& size_distribution,
  bool& use_cross_sections,
  bool& compute_phase_function, 
  bool& use_legendre_series, 
  std::vector<double>& phase_function_angles, 
  size_t& nb_legendre_terms,
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

  //the particle size is either a single radius or the description of a size distribution
  //additional comments and empty lines in front of it are skipped
  while (std::getline(file, line))
  {
    const size_t first = line.find_first_not_of(" \t\r");

    if (first != std::string::npos && line[first] != '#') break;
  }

  if (!readParticleSize(line, radius, use_size_distribution, size_distribution))
    return false;

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


      compute_phase_function = readPhaseFunctionParameter(
        file, 
        use_legendre_series, 
        phase_function_angles, 
        nb_legendre_terms, 
        phase_function_output_file);

      if (!compute_phase_function)
        std::cout << "Error reading the phase function/Legendre series parameters!\n";

    }
    else
      compute_phase_function = false;
  }


  file.close();


  std::cout << "read from parameter file: \n";
  std::cout << "refractive index: " << refractive_index_file << "\n";
  if (use_size_distribution)
    std::cout << "particle size distribution: " << size_distribution.radius.size() << " radii between "
              << size_distribution.radius.front() << " and " << size_distribution.radius.back() << " micron\n";
  else
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




bool readRefractiveIndex(
  const std::string file_name, 
  std::vector<double> &wavelengths, 
  std::vector< std::complex<double> > &refractive_index)
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


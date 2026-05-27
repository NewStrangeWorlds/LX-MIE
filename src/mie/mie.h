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



#ifndef MIE_H_
#define MIE_H_

#include <vector>
#include <complex>


namespace lxmie
{

struct MieResult {
  double q_ext;
  double q_sca;
  double q_abs;
  double asymmetry_parameter;
};


MieResult Mie(
  std::complex<double> refractive_index, 
  double size_parameter);

MieResult Mie(
  std::complex<double> refractive_index, 
  double size_parameter,
  std::size_t nb_legendre_moments, 
  std::vector<double>& legendre_moments);

MieResult Mie(
  std::complex<double> refractive_index, 
  double size_parameter,
  const std::vector<double>& angular_grid,
  std::vector<std::complex<double>>& s1, 
  std::vector<std::complex<double>>& s2);

double phaseFunction(
  std::complex<double> s1, 
  std::complex<double> s2,
  double size_parameter, 
  double q_sca);

std::vector<double> phaseFunction(
  const std::vector<std::complex<double>>& s1,
  const std::vector<std::complex<double>>& s2,
  double size_parameter, 
  double q_sca);

}

#endif

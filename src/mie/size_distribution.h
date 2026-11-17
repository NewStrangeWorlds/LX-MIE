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



#ifndef SIZE_DISTRIBUTION_H_
#define SIZE_DISTRIBUTION_H_

#include <vector>
#include <complex>
#include <cstddef>


namespace lxmie
{

//default number of radius grid points used for the analytical size distributions
inline constexpr std::size_t default_nb_distribution_points = 500;


//A particle size distribution, discretised into a set of quadrature nodes
//  radius: particle radii
//  weight: number of particles represented by each node, i.e. n(r_i) * dr_i
//          for a continuous distribution n(r) = dN/dr
//The weights do not need to be normalised, all size-integrated results are
//normalised to the total number of particles (sum of all weights).
//The struct can be filled directly, e.g. with the size bins of a cloud model,
//or created with one of the functions below.
struct SizeDistribution {
  std::vector<double> radius;
  std::vector<double> weight;
};


//Size-integrated optical properties
//Cross sections are mean values per particle, in units of the squared
//length unit that is used for the wavelength and the particle radii
struct MieDistributionResult {
  double c_ext;                    //mean extinction cross section
  double c_sca;                    //mean scattering cross section
  double c_abs;                    //mean absorption cross section
  double q_ext;                    //effective efficiencies, i.e. the mean cross
  double q_sca;                    //sections divided by the mean geometric
  double q_abs;                    //cross section
  double asymmetry_parameter;      //scattering cross section-weighted mean
  double geometric_cross_section;  //mean geometric cross section <pi r^2>
  double effective_radius;         //<r^3> / <r^2>
};


//Log-normal distribution
//  n(r) ~ 1/r * exp( - ln^2(r/median_radius) / (2 ln^2(sigma_g)) )
//with the geometric standard deviation sigma_g > 1
//The radius grid is spaced equidistantly in ln(r). If r_min or r_max are
//zero, the corresponding integration limit is chosen automatically.
SizeDistribution logNormalDistribution(
  double median_radius,
  double sigma_g,
  std::size_t nb_points = default_nb_distribution_points,
  double r_min = 0.0,
  double r_max = 0.0);

//Gamma distribution following Hansen & Travis (1974)
//  n(r) ~ r^((1-3b)/b) * exp( - r/(ab) )
//with the effective radius a and the effective variance b (0 < b < 0.5)
//The radius grid is spaced equidistantly in ln(r). If r_min or r_max are
//zero, the corresponding integration limit is chosen automatically.
SizeDistribution gammaDistribution(
  double effective_radius,
  double effective_variance,
  std::size_t nb_points = default_nb_distribution_points,
  double r_min = 0.0,
  double r_max = 0.0);

//Power-law distribution n(r) ~ r^(-exponent) between r_min and r_max
//The radius grid is spaced equidistantly in ln(r).
SizeDistribution powerLawDistribution(
  double r_min,
  double r_max,
  double exponent,
  std::size_t nb_points = default_nb_distribution_points);

//Tabulated distribution n(r) = dN/dr, given at a set of strictly increasing radii
//The weights are obtained from the trapezoidal rule.
SizeDistribution tabulatedDistribution(
  const std::vector<double>& radius,
  const std::vector<double>& dn_dr);


//Mie calculations integrated over a particle size distribution
//The wavelength has to be given in the same units as the particle radii
//The sign of the imaginary part of the refractive index does not matter
//If compiled with OpenMP, the loop over the particle sizes is run in parallel,
//unless the function is itself called from within a parallel region.
MieDistributionResult MieDistribution(
  std::complex<double> refractive_index,
  double wavelength,
  const SizeDistribution& distribution);

//Additionally returns the Legendre series of the size-averaged phase function
MieDistributionResult MieDistribution(
  std::complex<double> refractive_index,
  double wavelength,
  const SizeDistribution& distribution,
  std::size_t nb_legendre_moments,
  std::vector<double>& legendre_moments);

//Additionally returns the size-averaged phase function
//on a grid of cosines of the scattering angle
MieDistributionResult MieDistribution(
  std::complex<double> refractive_index,
  double wavelength,
  const SizeDistribution& distribution,
  const std::vector<double>& angular_grid,
  std::vector<double>& phase_function);

}

#endif

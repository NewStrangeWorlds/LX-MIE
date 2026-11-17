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



#include "size_distribution.h"
#include "mie.h"

#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <stdexcept>
#include <atomic>
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace lxmie
{

inline constexpr double pi = 3.14159265358979323846;



//Creates a grid with nodes equidistantly spaced in ln(r) and the corresponding
//quadrature weights of the trapezoidal rule, i.e. dr = r * dln(r)
static SizeDistribution logarithmicGrid(
  const double r_min,
  const double r_max,
  const std::size_t nb_points)
{
  if (nb_points < 2)
    throw std::invalid_argument("LX-MIE: a size distribution requires at least two grid points");

  if (!(r_min > 0.0) || !(r_max > r_min) || !std::isfinite(r_max))
    throw std::invalid_argument("LX-MIE: the size distribution requires 0 < r_min < r_max");

  SizeDistribution distribution;
  distribution.radius.resize(nb_points);
  distribution.weight.resize(nb_points);

  const double log_min = std::log(r_min);
  const double step = (std::log(r_max) - log_min) / double(nb_points - 1);

  for (std::size_t i = 0; i < nb_points; ++i)
  {
    distribution.radius[i] = std::exp(log_min + double(i) * step);
    distribution.weight[i] = distribution.radius[i] * step;
  }

  distribution.radius.back() = r_max;

  distribution.weight.front() *= 0.5;
  distribution.weight.back()  *= 0.5;

  return distribution;
}




//The automatic integration limits are chosen such that they cover the
//distribution itself (required for the normalisation), as well as the
//distribution weighted with r^6. The latter is the steepest possible size
//dependence of a cross section (scattering in the Rayleigh limit).
SizeDistribution logNormalDistribution(
  const double median_radius,
  const double sigma_g,
  const std::size_t nb_points,
  const double r_min,
  const double r_max)
{
  if (!(median_radius > 0.0) || !std::isfinite(median_radius))
    throw std::invalid_argument("LX-MIE: log-normal distribution requires a median radius > 0");

  if (!(sigma_g > 1.0) || !std::isfinite(sigma_g))
    throw std::invalid_argument("LX-MIE: log-normal distribution requires a geometric standard deviation > 1");

  const double s = std::log(sigma_g);

  //r^k n(r) is a Gaussian in ln(r) with a width of s, centred at ln(median_radius) + k s^2
  const double nb_sigma = 6.0;
  const double lower = (r_min > 0.0) ? r_min : median_radius * std::exp(-nb_sigma * s);
  const double upper = (r_max > 0.0) ? r_max : median_radius * std::exp((6.0 * s + nb_sigma) * s);

  SizeDistribution distribution = logarithmicGrid(lower, upper, nb_points);

  for (std::size_t i = 0; i < nb_points; ++i)
  {
    const double z = std::log(distribution.radius[i] / median_radius) / s;

    distribution.weight[i] *= std::exp(-0.5 * z * z) / distribution.radius[i];
  }

  return distribution;
}




SizeDistribution gammaDistribution(
  const double effective_radius,
  const double effective_variance,
  const std::size_t nb_points,
  const double r_min,
  const double r_max)
{
  if (!(effective_radius > 0.0) || !std::isfinite(effective_radius))
    throw std::invalid_argument("LX-MIE: gamma distribution requires an effective radius > 0");

  if (!(effective_variance > 0.0) || !(effective_variance < 0.5))
    throw std::invalid_argument("LX-MIE: gamma distribution requires an effective variance between 0 and 0.5");

  //with t = r/(ab), the distribution is given by n(t) ~ t^q * exp(-t)
  const double scale = effective_radius * effective_variance;
  const double q = (1.0 - 3.0 * effective_variance) / effective_variance;

  //lower limit: the neglected fraction of the total particle number is about 1e-8
  double log_t_min = (std::log(1e-8) + std::lgamma(q + 2.0)) / (q + 1.0);
  if (log_t_min < std::log(1e-8)) log_t_min = std::log(1e-8);

  //upper limit: far in the exponential tail of t^(q+6) * exp(-t)
  const double p = q + 6.0;
  const double t_max = p + 10.0 * std::sqrt(p + 1.0) + 20.0;

  const double lower = (r_min > 0.0) ? r_min : scale * std::exp(log_t_min);
  const double upper = (r_max > 0.0) ? r_max : scale * t_max;

  SizeDistribution distribution = logarithmicGrid(lower, upper, nb_points);

  //the distribution is scaled with its maximum value to avoid overflows
  const double log_max = (q > 0.0) ? q * std::log(q) - q : 0.0;

  for (std::size_t i = 0; i < nb_points; ++i)
  {
    const double t = distribution.radius[i] / scale;

    distribution.weight[i] *= std::exp(q * std::log(t) - t - log_max);
  }

  return distribution;
}




SizeDistribution powerLawDistribution(
  const double r_min,
  const double r_max,
  const double exponent,
  const std::size_t nb_points)
{
  if (!std::isfinite(exponent))
    throw std::invalid_argument("LX-MIE: power-law distribution requires a finite exponent");

  SizeDistribution distribution = logarithmicGrid(r_min, r_max, nb_points);

  //scaled to n(r_min) = 1
  for (std::size_t i = 0; i < nb_points; ++i)
    distribution.weight[i] *= std::exp(-exponent * std::log(distribution.radius[i] / r_min));

  return distribution;
}




SizeDistribution tabulatedDistribution(
  const std::vector<double>& radius,
  const std::vector<double>& dn_dr)
{
  if (radius.size() != dn_dr.size())
    throw std::invalid_argument("LX-MIE: tabulated distribution requires the same number of radii and distribution values");

  if (radius.size() < 2)
    throw std::invalid_argument("LX-MIE: tabulated distribution requires at least two radii");

  const std::size_t nb_points = radius.size();

  for (std::size_t i = 1; i < nb_points; ++i)
    if (!(radius[i] > radius[i-1]))
      throw std::invalid_argument("LX-MIE: the radii of a tabulated distribution have to be strictly increasing");

  SizeDistribution distribution;
  distribution.radius = radius;
  distribution.weight.assign(nb_points, 0.0);

  //trapezoidal rule
  for (std::size_t i = 0; i < nb_points - 1; ++i)
  {
    const double half_step = 0.5 * (radius[i+1] - radius[i]);

    distribution.weight[i]   += half_step * dn_dr[i];
    distribution.weight[i+1] += half_step * dn_dr[i+1];
  }

  return distribution;
}




static void checkDistribution(const SizeDistribution& distribution)
{
  if (distribution.radius.size() != distribution.weight.size())
    throw std::invalid_argument("LX-MIE: size distribution requires the same number of radii and weights");

  if (distribution.radius.empty())
    throw std::invalid_argument("LX-MIE: size distribution is empty");

  double weight_sum = 0.0;

  for (std::size_t i = 0; i < distribution.radius.size(); ++i)
  {
    if (!(distribution.radius[i] > 0.0) || !std::isfinite(distribution.radius[i]))
      throw std::invalid_argument("LX-MIE: size distribution radii have to be positive");

    if (!(distribution.weight[i] >= 0.0) || !std::isfinite(distribution.weight[i]))
      throw std::invalid_argument("LX-MIE: size distribution weights have to be non-negative");

    weight_sum += distribution.weight[i];
  }

  if (!(weight_sum > 0.0) || !std::isfinite(weight_sum))
    throw std::invalid_argument("LX-MIE: the sum of the size distribution weights has to be positive");
}




enum class AngularOutput { none, legendre_series, phase_function };


//Performs the Mie calculations for all particle sizes and integrates the results
//over the size distribution. Depending on angular_output, angular_result either
//contains the Legendre series of the phase function or the phase function at
//the angles given by angular_grid. Both are averaged with the scattering cross
//sections as weights.
static MieDistributionResult integrateDistribution(
  const std::complex<double> refractive_index,
  const double wavelength,
  const SizeDistribution& distribution,
  const AngularOutput angular_output,
  const std::size_t nb_legendre_moments,
  const std::vector<double>* angular_grid,
  std::vector<double>* angular_result)
{
  checkDistribution(distribution);

  if (!(wavelength > 0.0) || !std::isfinite(wavelength))
    throw std::invalid_argument("LX-MIE: the wavelength has to be positive");

  const std::size_t nb_sizes = distribution.radius.size();

  std::size_t nb_angular = 0;
  if (angular_output == AngularOutput::legendre_series) nb_angular = nb_legendre_moments;
  if (angular_output == AngularOutput::phase_function)  nb_angular = angular_grid->size();

  //cross sections of the individual particle sizes
  std::vector<double> c_ext(nb_sizes, 0.0);
  std::vector<double> c_sca(nb_sizes, 0.0);
  std::vector<double> asymmetry(nb_sizes, 0.0);

  //angular quantities of the individual sizes, multiplied by their scattering cross sections
  std::vector<std::vector<double>> angular(nb_angular > 0 ? nb_sizes : 0);

  std::atomic<bool> error_flag{false};
  std::string error_message;

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) if(!omp_in_parallel())
  #endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(nb_sizes); ++j)
  {
    const std::size_t i = static_cast<std::size_t>(j);

    if (error_flag || distribution.weight[i] == 0.0) continue;

    try
    {
      const double radius = distribution.radius[i];
      const double size_parameter = 2.0 * pi * radius / wavelength;
      const double area = pi * radius * radius;

      MieResult r{};

      if (angular_output == AngularOutput::legendre_series)
      {
        angular[i].assign(nb_angular, 0.0);
        r = Mie(refractive_index, size_parameter, nb_angular, angular[i]);

        //the moments returned by Mie are normalised by the scattering efficiency
        for (auto& moment : angular[i])
          moment = (r.q_sca > 0.0) ? moment * r.q_sca * area : 0.0;
      }
      else if (angular_output == AngularOutput::phase_function)
      {
        thread_local std::vector<std::complex<double>> s1, s2;
        s1.resize(nb_angular);
        s2.resize(nb_angular);

        r = Mie(refractive_index, size_parameter, *angular_grid, s1, s2);

        //phase function times scattering cross section, see Eq. 7
        const double factor = 2.0 * area / (size_parameter * size_parameter);

        angular[i].resize(nb_angular);
        for (std::size_t k = 0; k < nb_angular; ++k)
          angular[i][k] = factor * (std::norm(s1[k]) + std::norm(s2[k]));
      }
      else
        r = Mie(refractive_index, size_parameter);

      c_ext[i] = r.q_ext * area;
      c_sca[i] = r.q_sca * area;
      asymmetry[i] = r.asymmetry_parameter;
    }
    catch (const std::exception& e)
    {
      #ifdef _OPENMP
      #pragma omp critical(lxmie_distribution_error)
      #endif
      { error_flag = true; error_message = e.what(); }
    }
  }

  if (error_flag) throw std::runtime_error(error_message);


  //integration over the size distribution
  double sum_weight = 0.0;
  double sum_r2 = 0.0;
  double sum_r3 = 0.0;
  double sum_ext = 0.0;
  double sum_sca = 0.0;
  double sum_asymmetry = 0.0;

  std::vector<double> sum_angular(nb_angular, 0.0);

  for (std::size_t i = 0; i < nb_sizes; ++i)
  {
    const double weight = distribution.weight[i];
    const double radius = distribution.radius[i];

    if (weight == 0.0) continue;

    sum_weight    += weight;
    sum_r2        += weight * radius * radius;
    sum_r3        += weight * radius * radius * radius;
    sum_ext       += weight * c_ext[i];
    sum_sca       += weight * c_sca[i];
    sum_asymmetry += weight * c_sca[i] * asymmetry[i];

    for (std::size_t k = 0; k < nb_angular; ++k)
      sum_angular[k] += weight * angular[i][k];
  }

  MieDistributionResult result{};

  result.geometric_cross_section = pi * sum_r2 / sum_weight;
  result.effective_radius = sum_r3 / sum_r2;

  result.c_ext = sum_ext / sum_weight;
  result.c_sca = sum_sca / sum_weight;
  result.c_abs = result.c_ext - result.c_sca;

  result.q_ext = result.c_ext / result.geometric_cross_section;
  result.q_sca = result.c_sca / result.geometric_cross_section;
  result.q_abs = result.c_abs / result.geometric_cross_section;

  result.asymmetry_parameter = (sum_sca > 0.0) ? sum_asymmetry / sum_sca : 0.0;

  if (angular_result != nullptr)
  {
    angular_result->assign(nb_angular, 0.0);

    if (sum_sca > 0.0)
      for (std::size_t k = 0; k < nb_angular; ++k)
        (*angular_result)[k] = sum_angular[k] / sum_sca;
  }

  return result;
}




MieDistributionResult MieDistribution(
  const std::complex<double> refractive_index,
  const double wavelength,
  const SizeDistribution& distribution)
{
  return integrateDistribution(
    refractive_index, wavelength, distribution,
    AngularOutput::none, 0, nullptr, nullptr);
}




MieDistributionResult MieDistribution(
  const std::complex<double> refractive_index,
  const double wavelength,
  const SizeDistribution& distribution,
  const std::size_t nb_legendre_moments,
  std::vector<double>& legendre_moments)
{
  return integrateDistribution(
    refractive_index, wavelength, distribution,
    AngularOutput::legendre_series, nb_legendre_moments, nullptr, &legendre_moments);
}




MieDistributionResult MieDistribution(
  const std::complex<double> refractive_index,
  const double wavelength,
  const SizeDistribution& distribution,
  const std::vector<double>& angular_grid,
  std::vector<double>& phase_function)
{
  return integrateDistribution(
    refractive_index, wavelength, distribution,
    AngularOutput::phase_function, 0, &angular_grid, &phase_function);
}




}

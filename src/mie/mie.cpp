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



#include "mie.h"

#include <iostream>
#include <complex>
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>


namespace lxmie
{

inline constexpr std::size_t continued_fraction_max_terms = 10000000;
inline constexpr double continued_fraction_epsilon = 1e-10;



//Calculates single terms a_n used in the evaluation of the continued fractions by Lentz (1976)
//See Eq. 9 in Lentz (1976)
template<typename T>
T anCoeff(const std::size_t n, const double nu, const T& z)
{
  return (n % 2 == 0 ? -1.0 : 1.0) * 2.0 * (nu + n - 1) / z;
}




//Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
//Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
//Returns A_N
//This is the version for a complex A_N
std::complex<double> startingANcontinuedFractions(const std::size_t nb_mie_terms, const std::complex<double> mx)
{
  const double nu = nb_mie_terms + 0.5;

  std::complex<double> function_numerator(1, 1);
  std::complex<double> function_denominator(1, 1);

  //n=1
  std::complex<double> a_numerator = anCoeff(std::size_t{1}, nu, mx);
  std::complex<double> a_denominator(1.0, 0.0);

  function_numerator *= a_numerator;
  function_denominator *= a_denominator;

  //n=2
  a_numerator   = anCoeff(std::size_t{2}, nu, mx) + 1.0 / a_numerator;
  a_denominator = anCoeff(std::size_t{2}, nu, mx);

  function_numerator *= a_numerator;
  function_denominator *= a_denominator;

  for (std::size_t i = 3; i < continued_fraction_max_terms; ++i)
  {
    const std::complex<double> a_i = anCoeff(i, nu, mx);

    a_numerator   = a_i + 1.0 / a_numerator;
    a_denominator = a_i + 1.0 / a_denominator;

    function_numerator   *= a_numerator;
    function_denominator *= a_denominator;

    if (std::abs(a_numerator - a_denominator) / std::abs(a_numerator) < continued_fraction_epsilon)
      return function_numerator / function_denominator - double(nb_mie_terms) / mx;
  }

  throw std::runtime_error("LX-MIE: complex continued fraction did not converge");
}




//Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
//Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
//Returns A_N
//This is the version for a real A_N
double startingANcontinuedFractionsReal(const std::size_t nb_mie_terms, const double x)
{
  const double nu = nb_mie_terms + 0.5;

  double function_numerator = 1.0;
  double function_denominator = 1.0;

  //n=1
  double a_numerator   = anCoeff(std::size_t{1}, nu, x);
  double a_denominator = 1.0;

  function_numerator   *= a_numerator;
  function_denominator *= a_denominator;

  //n=2
  a_numerator   = anCoeff(std::size_t{2}, nu, x) + 1.0 / a_numerator;
  a_denominator = anCoeff(std::size_t{2}, nu, x);

  function_numerator   *= a_numerator;
  function_denominator *= a_denominator;

  for (std::size_t i = 3; i < continued_fraction_max_terms; ++i)
  {
    const double a_i = anCoeff(i, nu, x);

    a_numerator   = a_i + 1.0 / a_numerator;
    a_denominator = a_i + 1.0 / a_denominator;

    function_numerator   *= a_numerator;
    function_denominator *= a_denominator;

    if (std::abs(a_numerator - a_denominator) / std::abs(a_numerator) < continued_fraction_epsilon)
      return function_numerator / function_denominator - double(nb_mie_terms) / x;
  }

  throw std::runtime_error("LX-MIE: real continued fraction did not converge");
}




//Calculates the Mie coefficients a and b that are used for the construction of the Mie series, see Eq. 17
//The required coefficients A_n are calculated via downward recursion, B_n and C_n by upward recursion
//a and b are evaluated up to a number of "nb_mie_terms"
void calcMieCoefficients(const std::size_t nb_mie_terms, const std::complex<double> refractive_index, const double size_parameter,
                         std::vector<std::complex<double>>& mie_coeff_a, std::vector<std::complex<double>>& mie_coeff_b)
{
  const double x = size_parameter;
  const std::complex<double> m = refractive_index;
  const std::complex<double> mx = m * x;

  //First, calculate A_n via backward recursion
  //A_n(mx) is complex; A_n(x) is real — stored separately
  std::vector<double> A_n_real(nb_mie_terms + 1);

  mie_coeff_a.back() = startingANcontinuedFractions(nb_mie_terms, mx);
  A_n_real.back()    = startingANcontinuedFractionsReal(nb_mie_terms, x);

  //backward recursion
  for (std::size_t n = nb_mie_terms; n > 1; --n)
  {
    mie_coeff_a[n-1] = double(n) / mx - 1.0 / (double(n) / mx + mie_coeff_a[n]);
    A_n_real[n-1]    = double(n) / x  - 1.0 / (double(n) / x  + A_n_real[n]);
  }


  //Now do a forward recursion to calculate B_n, C_n, and the Mie coefficients a_n and b_n
  std::complex<double> C_n(0.0, 0.0);
  std::complex<double> D_n(0.0, -1.0);

  //n=1
  C_n = std::complex<double>(1.0, (std::cos(x) + x * std::sin(x)) / (std::sin(x) - x * std::cos(x)));
  C_n = 1.0 / C_n;
  D_n = -1.0 / x + 1.0 / (1.0 / x - D_n);

  std::complex<double> A_n = mie_coeff_a[1];
  double A_nr = A_n_real[1];

  mie_coeff_a[1] = C_n * (A_n / m - A_nr) / (A_n / m - D_n);
  mie_coeff_b[1] = C_n * (A_n * m - A_nr) / (A_n * m - D_n);

  //n > 1
  for (std::size_t n = 2; n < mie_coeff_a.size(); ++n)
  {
    A_n  = mie_coeff_a[n];
    A_nr = A_n_real[n];

    D_n = -double(n) / x + 1.0 / (double(n) / x - D_n);
    C_n = C_n * (D_n + double(n) / x) / (A_nr + double(n) / x);

    mie_coeff_a[n] = C_n * (A_n / m - A_nr) / (A_n / m - D_n);
    mie_coeff_b[n] = C_n * (A_n * m - A_nr) / (A_n * m - D_n);
  }
}




//Calculates the Mie efficiencies, see Eq. 1
//The absorption efficiency is calculated as the difference of the extinction and scattering efficiencies
void calcMieEfficiencies(const double size_parameter,
                         const std::vector<std::complex<double>>& mie_coeff_a,
                         const std::vector<std::complex<double>>& mie_coeff_b,
                         double& q_ext, double& q_sca, double& q_abs)
{
  q_ext = 0.0;
  q_sca = 0.0;

  for (std::size_t n = 1; n < mie_coeff_a.size(); ++n)
  {
    q_sca += (2.0*n + 1.0) * (std::norm(mie_coeff_a[n]) + std::norm(mie_coeff_b[n]));
    q_ext += (2.0*n + 1.0) * std::real(mie_coeff_a[n] + mie_coeff_b[n]);
  }

  q_sca *= 2.0 / (size_parameter * size_parameter);
  q_ext *= 2.0 / (size_parameter * size_parameter);

  q_abs = q_ext - q_sca;
}




//Calculate and return the asymmetry parameter
//See Bohren&Huffman, page 120, for details on the equation
double calcAsymmetryParameter(const double q_sca, const double size_parameter,
                              const std::vector<std::complex<double>>& mie_coeff_a,
                              const std::vector<std::complex<double>>& mie_coeff_b)
{
  if (q_sca == 0.0) return 0.0;

  double g = 0.0;

  for (std::size_t n = 1; n < mie_coeff_a.size() - 1; ++n)
    g += n * (n + 2.0) / (n + 1.0) * std::real(mie_coeff_a[n] * std::conj(mie_coeff_a[n+1]) + mie_coeff_b[n] * std::conj(mie_coeff_b[n+1]))
         + (2.0*n + 1.0) / (n * (n + 1.0)) * std::real(mie_coeff_a[n] * std::conj(mie_coeff_b[n]));

  return g * 4.0 / (size_parameter * size_parameter * q_sca);
}




//Calculate the angular eigenfunction pi and tau at a specified angle (Eq. 6)
//Instead of evaluating the Legendre polynomials directly, we employ upward recurrence relations
//We here follow the notation of Wiscombe (1979)
//Note: cos_theta is the cosine of the scattering angle, not the angle itself
void calcAngularFunctions(const double cos_theta, std::vector<double>& pi_n, std::vector<double>& tau_n)
{
  pi_n[0] = 0.0;
  pi_n[1] = 1.0;

  double s, t;

  for (std::size_t n = 1; n < pi_n.size() - 1; ++n)
  {
    s = cos_theta * pi_n[n];
    t = s - pi_n[n-1];

    pi_n[n+1] = s + (n + 1.0) / double(n) * t;
    tau_n[n]  = double(n) * t - pi_n[n-1];
  }

  const std::size_t n = pi_n.size() - 1;
  s = cos_theta * pi_n[n];
  t = s - pi_n[n-1];
  tau_n[n] = double(n) * t - pi_n[n-1];
}




//Calculate the scattering amplitudes S_1 and S_2 at a specified angle
//The Mie intensities in Eq. 5 are given by i_1 = (abs[S_1])^2, i_2 = (abs[S_2])^2
//The intensities/amplitudes can be used to calculate the scattering phase function (see Eq. 7)
//See Wiscombe (1979) for the recurrence relations used below
//Note: cos_theta is the cosine of the scattering angle, not the angle itself
void calcScatteringAmplitudes(const double cos_theta,
                              const std::vector<std::complex<double>>& mie_coeff_a,
                              const std::vector<std::complex<double>>& mie_coeff_b,
                              std::complex<double>& s_1, std::complex<double>& s_2)
{
  std::vector<double> pi_n(mie_coeff_a.size(), 0.0);
  std::vector<double> tau_n(mie_coeff_a.size(), 0.0);

  calcAngularFunctions(cos_theta, pi_n, tau_n);

  std::complex<double> s_plus(0.0, 0.0);
  std::complex<double> s_minus(0.0, 0.0);

  for (std::size_t n = 1; n < pi_n.size(); ++n)
  {
    const double factor = (2.0*n + 1.0) / (double(n) * (n + 1.0));
    s_plus  += factor * (mie_coeff_a[n] + mie_coeff_b[n]) * (pi_n[n] + tau_n[n]);
    s_minus += factor * (mie_coeff_a[n] - mie_coeff_b[n]) * (pi_n[n] - tau_n[n]);
  }

  s_1 = 0.5 * (s_plus + s_minus);
  s_2 = 0.5 * (s_plus - s_minus);
}




//Calculate the Mueller arrays D and C, required for the Legendre series of the phase function
//This follows the procedure and notation of Dave (1988)
void calcMuellerArrays(const std::vector<std::complex<double>>& mie_coeff_a,
                       const std::vector<std::complex<double>>& mie_coeff_b,
                       const std::size_t nb_mie_terms,
                       std::vector<std::complex<double>>& C, std::vector<std::complex<double>>& D)
{
  if (nb_mie_terms < 2)
    throw std::runtime_error("LX-MIE: nb_mie_terms must be >= 2 for Mueller array calculation");

  C[nb_mie_terms + 2] = std::complex<double>(0.0, 0.0);
  D[nb_mie_terms + 2] = std::complex<double>(0.0, 0.0);

  C[nb_mie_terms + 1] = (1.0 - 1.0 / double(nb_mie_terms + 1)) * mie_coeff_b[nb_mie_terms];
  D[nb_mie_terms + 1] = (1.0 - 1.0 / double(nb_mie_terms + 1)) * mie_coeff_a[nb_mie_terms];

  C[nb_mie_terms] = (1.0 / double(nb_mie_terms) + 1.0 / double(nb_mie_terms + 1)) * mie_coeff_a[nb_mie_terms]
                  + (1.0 - 1.0 / double(nb_mie_terms)) * mie_coeff_b[nb_mie_terms - 1];
  D[nb_mie_terms] = (1.0 / double(nb_mie_terms) + 1.0 / double(nb_mie_terms + 1)) * mie_coeff_b[nb_mie_terms]
                  + (1.0 - 1.0 / double(nb_mie_terms)) * mie_coeff_a[nb_mie_terms - 1];

  for (int k = static_cast<int>(nb_mie_terms) - 1; k >= 2; --k)
  {
    C[k] = C[k+2] - (1.0 + 1.0 / double(k+1)) * mie_coeff_b[k+1]
                  + (1.0 / double(k) + 1.0 / double(k+1)) * mie_coeff_a[k]
                  + (1.0 - 1.0 / double(k)) * mie_coeff_b[k-1];

    D[k] = D[k+2] - (1.0 + 1.0 / double(k+1)) * mie_coeff_a[k+1]
                  + (1.0 / double(k) + 1.0 / double(k+1)) * mie_coeff_b[k]
                  + (1.0 - 1.0 / double(k)) * mie_coeff_a[k-1];
  }

  C[1] = C[3] + 1.5 * (mie_coeff_a[1] - mie_coeff_b[2]);
  D[1] = D[3] + 1.5 * (mie_coeff_b[1] - mie_coeff_a[2]);

  for (std::size_t k = 1; k <= nb_mie_terms + 2; ++k)
  {
    C[k] = double(2*k - 1) * C[k];
    D[k] = double(2*k - 1) * D[k];
  }
}




//Calculate the moments of the Legendre series of the phase function
//It follows the procedure of Dave (1988)
void calcLegendreMoments(const std::vector<std::complex<double>>& mie_coeff_a,
                         const std::vector<std::complex<double>>& mie_coeff_b,
                         const std::size_t nb_max_moments, const std::size_t nb_mie_terms,
                         std::vector<double>& legendre_moments)
{
  std::vector<double> a_m(nb_mie_terms + 2);
  std::vector<double> b_i(nb_max_moments + 2);
  std::vector<double> b_i_delta(nb_max_moments + 2);

  std::vector<std::complex<double>> C(nb_mie_terms + 5);
  std::vector<std::complex<double>> D(nb_mie_terms + 5);

  calcMuellerArrays(mie_coeff_a, mie_coeff_b, nb_mie_terms, C, D);

  const int N = static_cast<int>(nb_mie_terms);

  for (std::size_t l = 0; l < nb_max_moments; ++l)
  {
    const bool even_moment = !(l & 1);
    const int ld2 = static_cast<int>(l) / 2;
    int small_delta;

    legendre_moments[l] = 0.0;

    if (l == 0)
    {
      small_delta = 1;

      for (int m = 0; m <= N; ++m)
        a_m[m] = 2.0 / (2.0*m + 1.0);

      b_i[0] = 1.0;
    }
    else if (even_moment)
    {
      small_delta = 1;

      for (int m = ld2; m <= N; ++m)
        a_m[m] *= (1.0 + 1.0 / (2.0*m - static_cast<int>(l) + 1));

      for (int i = 0; i <= ld2 - 1; ++i)
        b_i[i] *= (1.0 - 1.0 / (static_cast<int>(l) - 2.0*i));

      b_i[ld2] = (2.0 - 1.0 / static_cast<int>(l)) * b_i[ld2 - 1];
    }
    else
    {
      small_delta = 2;

      for (int m = ld2; m <= N; ++m)
        a_m[m] *= (1.0 - 1.0 / (2.0*m + static_cast<int>(l) + 2.0));

      for (int i = 0; i <= ld2; ++i)
        b_i[i] *= (1.0 - 1.0 / (static_cast<int>(l) + 2.0*i + 1.0));
    }

    const int mmax = N - small_delta + 1;
    const int imax = std::min(ld2, mmax - ld2);

    for (int i = 0; i <= imax; ++i)
      b_i_delta[i] = b_i[i];

    if (even_moment) b_i_delta[0] *= 0.5;

    for (int i = 0; i <= imax; ++i)
    {
      double sum = 0.0;

      for (int m = ld2; m <= mmax - i; ++m)
        sum += a_m[m] * (std::real(C[m-i+1] * std::conj(C[m+i+small_delta]))
                       + std::real(D[m-i+1] * std::conj(D[m+i+small_delta])));

      legendre_moments[l] += b_i_delta[i] * sum;
    }

    legendre_moments[l] *= 0.5;
  }
}




//Calculate the maximum number of terms in the Mie series
//See Eq. 22
std::size_t numberOfMieTerms(const double size_parameter)
{
  return std::size_t(size_parameter + 4.3 * std::pow(size_parameter, 1.0/3.0) + 2);
}




//Internal helper: allocate coefficients, run the core calculation, and return a MieResult.
static MieResult calcMie(const std::complex<double> m, const double x,
                          std::vector<std::complex<double>>& a,
                          std::vector<std::complex<double>>& b)
{
  const std::size_t N = numberOfMieTerms(x);
  a.assign(N + 1, {});
  b.assign(N + 1, {});

  calcMieCoefficients(N, m, x, a, b);

  MieResult r{};
  calcMieEfficiencies(x, a, b, r.q_ext, r.q_sca, r.q_abs);
  r.asymmetry_parameter = calcAsymmetryParameter(r.q_sca, x, a, b);

  return r;
}




MieResult Mie(const std::complex<double> refractive_index, const double size_parameter)
{
  std::vector<std::complex<double>> a, b;
  return calcMie(refractive_index, size_parameter, a, b);
}




MieResult Mie(const std::complex<double> refractive_index, const double size_parameter,
              const std::size_t nb_legendre_moments, std::vector<double>& legendre_moments)
{
  std::vector<std::complex<double>> a, b;
  MieResult r = calcMie(refractive_index, size_parameter, a, b);

  const std::size_t N = a.size() - 1;
  std::vector<double> moments(nb_legendre_moments, 0.0);
  calcLegendreMoments(a, b, nb_legendre_moments, N, moments);

  for (std::size_t i = 0; i < nb_legendre_moments; ++i)
    legendre_moments[i] = 4.0 * moments[i] / (size_parameter * size_parameter * r.q_sca);

  return r;
}




MieResult Mie(const std::complex<double> refractive_index, const double size_parameter,
              const std::vector<double>& angular_grid,
              std::vector<std::complex<double>>& s1, std::vector<std::complex<double>>& s2)
{
  std::vector<std::complex<double>> a, b;
  MieResult r = calcMie(refractive_index, size_parameter, a, b);

  for (std::size_t i = 0; i < angular_grid.size(); ++i)
    calcScatteringAmplitudes(angular_grid[i], a, b, s1[i], s2[i]);

  return r;
}




//Returns the value of the scattering phase function for given (angular-dependent) scattering amplitudes
//See Eq. 7
double phaseFunction(const std::complex<double> s1, const std::complex<double> s2,
                     const double size_parameter, const double q_sca)
{
  return 2.0 / (size_parameter * size_parameter) * (std::norm(s1) + std::norm(s2)) / q_sca;
}




//Returns the values of the scattering phase function for different angles with given (angular-dependent) scattering amplitudes
//See Eq. 7
std::vector<double> phaseFunction(const std::vector<std::complex<double>>& s1,
                                  const std::vector<std::complex<double>>& s2,
                                  const double size_parameter, const double q_sca)
{
  std::vector<double> phase_function(s1.size(), 0.0);

  for (std::size_t i = 0; i < phase_function.size(); ++i)
    phase_function[i] = phaseFunction(s1[i], s2[i], size_parameter, q_sca);

  return phase_function;
}




}

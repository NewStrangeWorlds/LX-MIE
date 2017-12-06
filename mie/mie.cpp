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
#include <fstream>
#include <complex>
#include <algorithm>
#include <limits>
#include <cmath>


namespace lxmie
{

const size_t continued_fraction_max_terms = 10000001;
const double continued_fraction_epsilon = 1e-10;



//Calculates single terms a_n used in the evaluation of the continued fractions by Lentz (1976)
//See Eq. 9 in Lentz (1976)
//This is the version for a complex a_n
std::complex<double> an(const size_t n, const double nu, const std::complex<double>& z)
{

  if (n % 2 == 0)
    return -1.0 * 2. * (nu + n -1) * 1./z;
  else
    return 2. * (nu + n -1) * 1./z;

}




//Calculates single terms a_n used in the evaluation of the continued fractions by Lentz (1976)
//See Eq. 9 in Lentz (1976)
//This is the version for a real-numbered a_n
double anReal(const size_t n, const double nu, const double z)
{

  if (n % 2 == 0)
    return -1.0 * 2. * (nu + n -1) * 1./z;
  else
    return 2. * (nu + n -1) * 1./z;

}




//Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
//Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
//Returns A_N
//This is the version for a complex A_N
std::complex<double> startingANcontinuedFractions(const size_t nb_mie_terms, const std::complex<double> mx)
{
  double nu = nb_mie_terms + 0.5;

  //starting values
  std::complex<double> function_numerator(1,1);
  std::complex<double> function_denominator(1,1);


  //n=1
  std::complex<double> a_numerator = an(1, nu, mx);
  std::complex<double> a_demoninator = 1.0;

  function_numerator *= a_numerator;
  function_denominator *= a_demoninator;


  //n=2
  a_numerator = an(2, nu, mx) + 1./a_numerator;
  a_demoninator = an(2, nu, mx);

  function_numerator *= a_numerator;
  function_denominator *= a_demoninator;


  for (size_t i = 3; i < continued_fraction_max_terms; ++i)
  {
    std::complex<double> a_i = an(i, nu, mx);

    a_numerator = a_i + 1./a_numerator;
    a_demoninator = a_i + 1./a_demoninator;


    function_numerator *= a_numerator;
    function_denominator *= a_demoninator;


    if ( fabs( (std::abs(a_numerator) - std::abs(a_demoninator) ) / std::abs(a_numerator) ) < continued_fraction_epsilon)
      break;
  }


  return function_numerator / function_denominator - (nb_mie_terms*1.)/mx;
}




//Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
//Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
//Returns A_N
//This is the version for a real A_N
double startingANcontinuedFractionsReal(const size_t nb_mie_terms, const double x)
{
  double nu = nb_mie_terms + 0.5;


  //starting values
  double function_numerator = 1;
  double function_denominator = 1;


  //n = 1
  double a_numerator = anReal(1, nu, x);
  double a_denominator = 1.0;

  function_numerator *= a_numerator;
  function_denominator *= a_denominator;


  //n = 2
  a_numerator = anReal(2, nu, x) + 1./a_numerator;
  a_denominator = anReal(2, nu, x);

  function_numerator *= a_numerator;
  function_denominator *= a_denominator;


  for (size_t i = 3; i < continued_fraction_max_terms; ++i)
  {
    double a_i = anReal(i, nu, x);

    a_numerator = a_i + 1./a_numerator;
    a_denominator = a_i + 1./a_denominator;

    function_numerator *= a_numerator;
    function_denominator *= a_denominator;

    if ( std::abs( (a_numerator - a_denominator ) / a_numerator) < continued_fraction_epsilon )
      break;
  }



  return function_numerator / function_denominator - nb_mie_terms/x;
}




//Calculates the Mie coefficients a and b that are used for the construction of the Mie series, see Eq. 17
//The required coefficients A_n are calculated via downward recursion, B_n and C_n by upward recursion
//a and b are evaluated up to a number of "nb_mie_terms"
void calcMieCoefficients(const size_t nb_mie_terms, const std::complex<double> refractive_index, const double size_parameter,
                         std::vector< std::complex<double> >& mie_coeff_a, std::vector< std::complex<double> >& mie_coeff_b)
{
  double x = size_parameter;
  std::complex<double> m = refractive_index;
  std::complex<double> mx = m * x;


  //First, calculate A_n via backward recursion
  //Note that we need A_N(mx), a complex number, and A_N(x), a real number (see Eq. 17)
  //We use the Mie a-coefficient to store the complex values and the Mie b-coefficients for the real numbers
  mie_coeff_a.back() = startingANcontinuedFractions(nb_mie_terms, mx);
  mie_coeff_b.back().real(startingANcontinuedFractionsReal(nb_mie_terms, x));


  //backward recursion
  for (size_t n=nb_mie_terms; n>1; --n)
  {
    mie_coeff_a[n-1] = (n*1.)/mx - 1./( (1.*n)/mx + mie_coeff_a[n]);
    mie_coeff_b[n-1].real((n*1.)/x - 1./( (1.*n)/x + mie_coeff_b[n].real() ));
  }



  //Now we do a forward recursion to calculate B_n, C_n, and the Mie coefficients a_n and b_n
  std::complex<double> C_n(0.0, 0.0);
  std::complex<double> D_n(0.0, -1.0);


  //n=1
  C_n = std::complex<double>(1., (std::cos(x) + x * std::sin(x)) / (std::sin(x) - x * std::cos(x)) );
  C_n = 1./C_n;
  D_n = (-1.)/x + 1.0/((1.)/x - D_n);

  std::complex<double>A_n = mie_coeff_a[1];
  double A_n_r = mie_coeff_b[1].real();


  mie_coeff_a[1] = C_n * (A_n / m - A_n_r) / (A_n / m - D_n);
  mie_coeff_b[1] = C_n * (A_n * m - A_n_r) / (A_n * m - D_n);


  //n > 1
  for (size_t n=2; n<mie_coeff_a.size(); ++n)
  {
    A_n = mie_coeff_a[n];
    A_n_r = mie_coeff_b[n].real();

    D_n = (-1.*n)/x + 1.0/((1.*n)/x - D_n);
    C_n = C_n * ( D_n + (1.*n)/x )/(A_n_r + (1.*n)/x );


    mie_coeff_a[n] = C_n * (A_n / m - A_n_r) / (A_n / m - D_n);
    mie_coeff_b[n] = C_n * (A_n * m - A_n_r) / (A_n * m - D_n);
  }

}




//Calculates the Mie efficiencies, see Eq. 1
//The absorption efficiency is calculated as the difference of the extinction and scattering efficiencies
void calcMieEfficiencies(const double size_parameter,
                         const std::vector< std::complex<double> >& mie_coeff_a, const std::vector< std::complex<double> >& mie_coeff_b,
                         double& q_ext, double& q_sca, double& q_abs)
{
  q_ext = 0;
  q_sca = 0;

  for (size_t n=1; n<mie_coeff_a.size(); ++n)
  {
    q_sca += (2.*n + 1.) * (std::abs(mie_coeff_a[n]) * std::abs(mie_coeff_a[n]) + std::abs(mie_coeff_b[n]) * std::abs(mie_coeff_b[n]));
    q_ext += (2.*n + 1.) * std::real(mie_coeff_a[n] + mie_coeff_b[n]);
  }

  q_sca *= 2./(size_parameter*size_parameter);
  q_ext *= 2./(size_parameter*size_parameter);

  q_abs = q_ext - q_sca;
}




//Calculate and return the asymmetry parameter
//See Bohren&Huffman, page 120, for details on the equation
double calcAsymmetryParameter(const double q_sca, const double size_parameter,
                              const std::vector< std::complex<double> >& mie_coeff_a, const std::vector< std::complex<double> >& mie_coeff_b)
{
  double g = 0;

  for (size_t n=1; n<mie_coeff_a.size()-1; ++n)
    g += n*(n + 2.0) / (n + 1.0) * std::real(mie_coeff_a[n] * std::conj(mie_coeff_a[n+1]) + mie_coeff_b[n] * std::conj(mie_coeff_b[n+1]))
         + (2.*n + 1.0) / (n * (n + 1.0)) * std::real(mie_coeff_a[n] * std::conj(mie_coeff_b[n]));

  return g * 4./(size_parameter * size_parameter * q_sca);
}




//Calculate the angular eigenfunction pi and tau at a specified angle (Eq. 6)
//Instead of evaluating the Legendre polynomials directly, we employ upward recurrence relations
//We here follow the notation of Wiscombe (1979)
void calcAngularFunctions(const double angle, std::vector<double>& pi_n, std::vector<double>& tau_n)
{
  pi_n[0] = 0;
  pi_n[1] = 1;

  double s,t;

  for (size_t n=1; n<pi_n.size()-1; ++n)
  {
    s = angle * pi_n[n];
    t = s - pi_n[n-1];

    pi_n[n+1] = s + (n + 1.0)/(n*1.0) * t;

    tau_n[n] = n * t - pi_n[n-1];
  }

  size_t n = pi_n.size()-1;

  s = angle * pi_n[n];
  t = s - pi_n[n-1];

  tau_n[n] = n * t - pi_n[n-1];

}




//Calculate the scattering amplitudes S_1 and S_2 at a specified angle
//The Mie intensities in Eq. 5 are given by i_1 = (abs[S_1])^2, i_2 = (abs[S_2])^2
//The intensities/amplitudes can be used to calculate the scattering phase function (see Eq. 7)
//See Wiscombe (1979) for the recurrence relations used below
void calcScatteringAmplitudes(const double angle,
                              const std::vector< std::complex<double> >& mie_coeff_a, const std::vector< std::complex<double> >& mie_coeff_b,
                              std::complex<double>& s_1, std::complex<double>& s_2)
{
  std::vector<double> pi_n(mie_coeff_a.size(), 0);
  std::vector<double> tau_n(mie_coeff_a.size(), 0);


  calcAngularFunctions(angle, pi_n, tau_n);


  std::complex<double> s_plus(0,0);
  std::complex<double> s_minus(0,0);


  for (size_t n=1; n<pi_n.size(); ++n)
  {
    s_plus  += (2.*n + 1.0)/(n * (n + 1.0)) * (mie_coeff_a[n] + mie_coeff_b[n]) * (pi_n[n] + tau_n[n]);
    s_minus += (2.*n + 1.0)/(n * (n + 1.0)) * (mie_coeff_a[n] - mie_coeff_b[n]) * (pi_n[n] - tau_n[n]);
  }


  s_1 = 0.5 * (s_plus + s_minus);
  s_2 = 0.5 * (s_plus - s_minus);
}





//Calculate the Mueller arrays D and C, required for the Legendre series of the phase function
//This follows the procedure and notation of Dave (1988)
void calcMuellerArrays(const std::vector< std::complex<double> >& mie_coeff_a, const std::vector< std::complex<double> >& mie_coeff_b,
                       const size_t nb_mie_terms,
                       std::vector< std::complex<double> >& C, std::vector< std::complex<double> >& D)
{

  C[nb_mie_terms + 2] = std::complex<double>(0, 0);
  D[nb_mie_terms + 2] = std::complex<double>(0, 0);

  C[nb_mie_terms + 1] = (1. - 1./(nb_mie_terms+1)) * mie_coeff_b[nb_mie_terms];
  D[nb_mie_terms + 1] = (1. - 1./(nb_mie_terms+1)) * mie_coeff_a[nb_mie_terms];

  C[nb_mie_terms] = (1./nb_mie_terms + 1./(nb_mie_terms+1)) * mie_coeff_a[nb_mie_terms] + (1. - 1./nb_mie_terms) * mie_coeff_b[nb_mie_terms-1];
  D[nb_mie_terms] = (1./nb_mie_terms + 1./(nb_mie_terms+1)) * mie_coeff_b[nb_mie_terms] + (1. - 1./nb_mie_terms) * mie_coeff_a[nb_mie_terms-1];


  for (size_t k=nb_mie_terms-1; k>=2; k--)
  {
    C[k] = C[k+2] - (1. + 1./(k+1)) * mie_coeff_b[k+1] + (1./k + 1./(k+1)) * mie_coeff_a[k] + (1. - 1./k) * mie_coeff_b[k-1];

    D[k] = D[k+2] - (1. + 1./(k+1)) * mie_coeff_a[k+1] + (1./k + 1./(k+1)) * mie_coeff_b[k] + (1. - 1./k) * mie_coeff_a[k-1];
  }


  C[1] = C[3] + 1.5 * (mie_coeff_a[1] - mie_coeff_b[2]);
  D[1] = D[3] + 1.5 * (mie_coeff_b[1] - mie_coeff_a[2]);


  for (size_t k=1; k<=nb_mie_terms+2; ++k)
  {
    C[k] = (2.*k - 1) * C[k];
    D[k] = (2.*k - 1) * D[k];
  }

}





//Calculate the moments of the Legendre series of the phase function
//It follows the procedure of Dave (1988)
//Sorry for the messy code ^^
void calcLegendreMoments(const std::vector< std::complex<double> >& mie_coeff_a, const std::vector< std::complex<double> >& mie_coeff_b,
                         const size_t nb_max_moments, const size_t nb_mie_terms,
                         std::vector<double>& legendre_moments)
{

  std::vector<double> a_m(nb_mie_terms+2);
  std::vector<double> b_i(nb_max_moments+2);
  std::vector<double> b_i_delta(nb_max_moments+2);

  std::vector< std::complex<double> > C(nb_mie_terms+5);
  std::vector< std::complex<double> > D(nb_mie_terms+5);


  calcMuellerArrays(mie_coeff_a, mie_coeff_b, nb_mie_terms, C, D);


  for (size_t l=0; l<nb_max_moments; ++l)
  {
    bool even_moment = !(l & 1);
    size_t small_delta;
    size_t ld2 = l/2;

    legendre_moments[l] = 0;


    if (l==0)
    {
      small_delta = 1;

      for (size_t m=0; m<=nb_mie_terms; ++m)
        a_m[m] = 2.0 * 1./(2.*m + 1.);

      b_i[0] = 1.0;
    }
    else if (even_moment == true)
    {
      small_delta = 1;

      for (size_t m = ld2; m <= nb_mie_terms; ++m)
        a_m[m] = (1. + 1./(2.*m - l + 1)) * a_m[m];

      for (size_t i=0; i<=ld2-1; ++i)
        b_i[i] = (1. - 1./(l - 2.*i)) * b_i[i];

      b_i[ld2] = (2. - 1./l) * b_i[ld2-1];
    }
    else if (even_moment == false)
    {
      small_delta = 2;

      for (size_t m = ld2; m <= nb_mie_terms; ++m)
        a_m[m] = (1. - 1./(2.*m + l + 2.)) * a_m[m];

      for (size_t i=0; i<=ld2; ++i)
        b_i[i] = (1. - 1./(l + 2.*i + 1)) * b_i[i];
    }

    int mmax = nb_mie_terms - small_delta;
    mmax++;

    int imax = std::min(ld2, mmax - ld2);

    for (int i=0; i<=imax; ++i)
      b_i_delta[i] = b_i[i];

    if (even_moment == true) b_i_delta[0] = 0.5 * b_i_delta[0];


    for (int i=0; i<=imax; ++i)
    {
      double sum = 0;

      for (int m=ld2; m <= mmax - i; m++)
        sum += a_m[m] * (std::real( C[m-i+1] * std::conj(C[m+i+small_delta]) ) + std::real( D[m-i+1] * std::conj(D[m+i+small_delta]) ));


      legendre_moments[l] += b_i_delta[i] * sum;
    }


    legendre_moments[l] *= 0.5;
  }

}




//Calculate the maximum number of terms in the Mie series
//See Eq. 22
size_t numberOfMieTerms(const double size_parameter)
{

  return size_t ( size_parameter + 4.3 * std::pow(size_parameter,1./3.) + 2);

}





int Mie(const std::complex<double> refractive_index, const double size_parameter, const size_t nb_legendre_moments,
        double& q_ext, double& q_sca, double& q_abs, double& asymmetry_parameter, std::vector<double>& legendre_moments)
{
  size_t nb_mie_terms =  numberOfMieTerms(size_parameter);


  std::vector< std::complex<double> > mie_coeff_a(nb_mie_terms+1);
  std::vector< std::complex<double> > mie_coeff_b(nb_mie_terms+1);


  calcMieCoefficients(nb_mie_terms, refractive_index, size_parameter, mie_coeff_a, mie_coeff_b);
  calcMieEfficiencies(size_parameter, mie_coeff_a, mie_coeff_b, q_ext, q_sca, q_abs);

  asymmetry_parameter = calcAsymmetryParameter(q_sca, size_parameter, mie_coeff_a, mie_coeff_b);


  std::vector<double> moments(nb_legendre_moments, 0);

  calcLegendreMoments(mie_coeff_a, mie_coeff_b, nb_legendre_moments, nb_mie_terms, moments);


  for (unsigned int i=0; i<moments.size(); i++)
    legendre_moments[i] = 4.*moments[i]/(size_parameter*size_parameter * q_sca);


  return 0;
}




int Mie(const std::complex<double> refractive_index, const double size_parameter, const std::vector<double>& angular_grid,
        double& q_ext, double& q_sca, double& q_abs, double& asymmetry_parameter,
        std::vector< std::complex<double> >& s1, std::vector< std::complex<double> >& s2)
{
  size_t nb_mie_terms =  numberOfMieTerms(size_parameter);


  std::vector< std::complex<double> > mie_coeff_a(nb_mie_terms+1);
  std::vector< std::complex<double> > mie_coeff_b(nb_mie_terms+1);


  calcMieCoefficients(nb_mie_terms, refractive_index, size_parameter, mie_coeff_a, mie_coeff_b);
  calcMieEfficiencies(size_parameter, mie_coeff_a, mie_coeff_b, q_ext, q_sca, q_abs);

  asymmetry_parameter = calcAsymmetryParameter(q_sca, size_parameter, mie_coeff_a, mie_coeff_b);


  for (size_t angle=0; angle<angular_grid.size(); ++angle)
  {
    std::complex<double> s_1;
    std::complex<double> s_2;

    calcScatteringAmplitudes(angular_grid[angle], mie_coeff_a, mie_coeff_b, s_1,s_2);

    s1[angle] = s_1;
    s2[angle] = s_2;
  }


  return 0;
}




int Mie(const std::complex<double> refractive_index, const double size_parameter,
        double& q_ext, double& q_sca, double& q_abs, double& asymmetry_parameter)
{
  size_t nb_mie_terms =  numberOfMieTerms(size_parameter);


  std::vector< std::complex<double> > mie_coeff_a(nb_mie_terms+1);
  std::vector< std::complex<double> > mie_coeff_b(nb_mie_terms+1);


  calcMieCoefficients(nb_mie_terms, refractive_index, size_parameter, mie_coeff_a, mie_coeff_b);
  calcMieEfficiencies(size_parameter, mie_coeff_a, mie_coeff_b, q_ext, q_sca, q_abs);

  asymmetry_parameter = calcAsymmetryParameter(q_sca, size_parameter, mie_coeff_a, mie_coeff_b);


  return 0;
}





//Returns the value of the scattering phase functions for given (angular-dependent) scattering amplitudes
//See Eq. 7
double phaseFunction(const std::complex<double> s1, const std::complex<double> s2, const double size_parameter, const double q_sca)
{

  return 2.0/(size_parameter*size_parameter) * (std::abs(s1)*std::abs(s1) + std::abs(s2) * std::abs(s2)) / q_sca;

}





//Returns the values of the scattering phase functions for different angles with given (angular-dependent) scattering amplitudes
//See Eq. 7
std::vector<double> phaseFunction(const std::vector< std::complex<double> >& s1, const std::vector< std::complex<double> >& s2,
                                  const double size_parameter, const double q_sca)
{
  std::vector<double> phase_function(s1.size(), 0);


  for (size_t i=0; i<phase_function.size(); ++i)
    phase_function[i] = phaseFunction(s1[i], s2[i], size_parameter, q_sca);


  return phase_function;
}




}



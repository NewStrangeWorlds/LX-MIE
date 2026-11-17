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

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include <atomic>
#include <cmath>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../mie/mie.h"
#include "../mie/size_distribution.h"

namespace py = pybind11;
using cd = std::complex<double>;


//LX-MIE internally uses a negative imaginary part of the refractive index,
//so the sign is enforced here to accept either convention (n + ik or n - ik)
inline cd refractiveIndex(const double n_real, const double n_imag)
{
  return cd(n_real, -std::abs(n_imag));
}


//Runs function(i) for i = 0 ... nb_elements-1 in an OpenMP-parallel loop with the GIL released
//Exceptions thrown within the loop are collected and re-thrown after the loop has finished
template <typename Function>
void parallelLoop(const std::size_t nb_elements, Function&& function)
{
  std::atomic<bool> error_flag{false};
  std::string error_message;
  {
    py::gil_scoped_release release;
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(nb_elements); ++i)
    {
      if (error_flag) continue;
      try
      {
        function(static_cast<std::size_t>(i));
      }
      catch (const std::exception& e)
      {
        #pragma omp critical
        { error_flag = true; error_message = e.what(); }
      }
    }
  }
  if (error_flag) throw std::runtime_error(error_message);
}


//The batch versions of the size distribution functions require lists of equal lengths
inline void checkListSizes(
  const std::vector<double>& n_real,
  const std::vector<double>& n_imag,
  const std::vector<double>& wavelengths)
{
  if (n_real.size() != n_imag.size() || n_real.size() != wavelengths.size())
    throw std::invalid_argument("LX-MIE: n_real, n_imag, and wavelength need to have the same length");
}


PYBIND11_MODULE(lxmie, m)
{
  m.doc() = "LX-MIE: Mie scattering for arbitrarily large size parameters";

  // -------------------------------------------------------------------------
  // MieResult
  // -------------------------------------------------------------------------

  py::class_<lxmie::MieResult>(m, "MieResult")
    .def_readonly("q_ext",               &lxmie::MieResult::q_ext)
    .def_readonly("q_sca",               &lxmie::MieResult::q_sca)
    .def_readonly("q_abs",               &lxmie::MieResult::q_abs)
    .def_readonly("asymmetry_parameter", &lxmie::MieResult::asymmetry_parameter)
    .def("__repr__", [](const lxmie::MieResult& r) {
      return "MieResult(q_ext=" + std::to_string(r.q_ext)
           + ", q_sca=" + std::to_string(r.q_sca)
           + ", q_abs=" + std::to_string(r.q_abs)
           + ", asymmetry_parameter=" + std::to_string(r.asymmetry_parameter) + ")";
    });


  // -------------------------------------------------------------------------
  // mie — basic Mie calculation (efficiencies + asymmetry parameter)
  //   Overload 1: scalar arguments → single MieResult
  //   Overload 2: list arguments  → list[MieResult], computed in parallel
  // -------------------------------------------------------------------------

  m.def("mie",
    [](double n_real, double n_imag, double x) {
      return lxmie::Mie(refractiveIndex(n_real, n_imag), x);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"),
    "Compute Mie efficiencies and asymmetry parameter for a single particle.");

  m.def("mie",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& size_params) {
      const std::size_t N = n_real.size();
      std::vector<lxmie::MieResult> results(N);
      std::atomic<bool> error_flag{false};
      std::string error_message;
      {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i)
        {
          if (error_flag) continue;
          try
          {
            results[i] = lxmie::Mie(refractiveIndex(n_real[i], n_imag[i]), size_params[i]);
          }
          catch (const std::runtime_error& e)
          {
            #pragma omp critical
            { error_flag = true; error_message = e.what(); }
          }
        }
      }
      if (error_flag) throw std::runtime_error(error_message);
      return results;
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"),
    "Compute Mie efficiencies and asymmetry parameter for a list of particles (OpenMP-parallel).");


  // -------------------------------------------------------------------------
  // mie_legendre — Mie calculation with Legendre series of the phase function
  //   Overload 1: scalar → (MieResult, list[float])
  //   Overload 2: lists  → (list[MieResult], list[list[float]])
  // -------------------------------------------------------------------------

  m.def("mie_legendre",
    [](double n_real, double n_imag, double x, std::size_t nb_moments) {
      std::vector<double> moments(nb_moments, 0.0);
      auto r = lxmie::Mie(refractiveIndex(n_real, n_imag), x, nb_moments, moments);
      return py::make_tuple(r, moments);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"), py::arg("nb_moments"),
    "Compute Mie efficiencies plus Legendre series of the phase function for a single particle.");

  m.def("mie_legendre",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& size_params,
       std::size_t nb_moments) {
      const std::size_t N = n_real.size();
      std::vector<lxmie::MieResult> results(N);
      std::vector<std::vector<double>> all_moments(N);
      std::atomic<bool> error_flag{false};
      std::string error_message;
      {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i)
        {
          if (error_flag) continue;
          try
          {
            all_moments[i].assign(nb_moments, 0.0);
            results[i] = lxmie::Mie(
              refractiveIndex(n_real[i], n_imag[i]), size_params[i], nb_moments, all_moments[i]);
          }
          catch (const std::runtime_error& e)
          {
            #pragma omp critical
            { error_flag = true; error_message = e.what(); }
          }
        }
      }
      if (error_flag) throw std::runtime_error(error_message);
      return py::make_tuple(results, all_moments);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"), py::arg("nb_moments"),
    "Compute Mie efficiencies plus Legendre series of the phase function for a list of particles (OpenMP-parallel).");


  // -------------------------------------------------------------------------
  // mie_amplitudes — Mie calculation with scattering amplitudes at given angles
  //   Overload 1: scalar n/x → (MieResult, list[complex], list[complex])
  //   Overload 2: lists  n/x → (list[MieResult], list[list[complex]], list[list[complex]])
  // -------------------------------------------------------------------------

  m.def("mie_amplitudes",
    [](double n_real, double n_imag, double x, const std::vector<double>& angles) {
      std::vector<cd> s1(angles.size()), s2(angles.size());
      auto r = lxmie::Mie(refractiveIndex(n_real, n_imag), x, angles, s1, s2);
      return py::make_tuple(r, s1, s2);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"), py::arg("angles"),
    "Compute Mie efficiencies plus scattering amplitudes at the given angles for a single particle.");

  m.def("mie_amplitudes",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& size_params,
       const std::vector<double>& angles) {
      const std::size_t N = n_real.size();
      std::vector<lxmie::MieResult> results(N);
      std::vector<std::vector<cd>> all_s1(N), all_s2(N);
      std::atomic<bool> error_flag{false};
      std::string error_message;
      {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(N); ++i)
        {
          if (error_flag) continue;
          try
          {
            all_s1[i].resize(angles.size());
            all_s2[i].resize(angles.size());
            results[i] = lxmie::Mie(
              refractiveIndex(n_real[i], n_imag[i]), size_params[i], angles, all_s1[i], all_s2[i]);
          }
          catch (const std::runtime_error& e)
          {
            #pragma omp critical
            { error_flag = true; error_message = e.what(); }
          }
        }
      }
      if (error_flag) throw std::runtime_error(error_message);
      return py::make_tuple(results, all_s1, all_s2);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("size_parameter"), py::arg("angles"),
    "Compute Mie efficiencies plus scattering amplitudes at the given angles for a list of particles (OpenMP-parallel).");


  // -------------------------------------------------------------------------
  // phase_function
  //   Overload 1: single s1/s2 scalars → float
  //   Overload 2: lists of s1/s2       → list[float]
  // -------------------------------------------------------------------------

  m.def("phase_function",
    py::overload_cast<cd, cd, double, double>(&lxmie::phaseFunction),
    py::arg("s1"), py::arg("s2"), py::arg("size_parameter"), py::arg("q_sca"),
    "Compute the phase function for a single scattering angle.");

  m.def("phase_function",
    py::overload_cast<const std::vector<cd>&, const std::vector<cd>&,
                      double, double>(&lxmie::phaseFunction),
    py::arg("s1"), py::arg("s2"), py::arg("size_parameter"), py::arg("q_sca"),
    "Compute the phase function for a list of scattering angles.");


  // -------------------------------------------------------------------------
  // SizeDistribution — particle size distribution as quadrature nodes/weights
  // -------------------------------------------------------------------------

  py::class_<lxmie::SizeDistribution>(m, "SizeDistribution")
    .def(py::init([](const std::vector<double>& radius, const std::vector<double>& weight) {
      if (radius.size() != weight.size())
        throw std::invalid_argument("LX-MIE: size distribution requires the same number of radii and weights");
      return lxmie::SizeDistribution{radius, weight};
    }),
    py::arg("radius"), py::arg("weight"),
    "Size distribution from a list of radii and the number of particles (weight) for each radius.")
    .def_readwrite("radius", &lxmie::SizeDistribution::radius)
    .def_readwrite("weight", &lxmie::SizeDistribution::weight)
    .def("__len__", [](const lxmie::SizeDistribution& d) { return d.radius.size(); })
    .def("__repr__", [](const lxmie::SizeDistribution& d) {
      if (d.radius.empty()) return std::string("SizeDistribution(empty)");
      return "SizeDistribution(" + std::to_string(d.radius.size()) + " radii, r_min="
           + std::to_string(d.radius.front()) + ", r_max=" + std::to_string(d.radius.back()) + ")";
    });

  m.def("lognormal_distribution", &lxmie::logNormalDistribution,
    py::arg("median_radius"), py::arg("sigma_g"),
    py::arg("nb_points") = lxmie::default_nb_distribution_points,
    py::arg("r_min") = 0.0, py::arg("r_max") = 0.0,
    "Log-normal size distribution with a median radius and a geometric standard deviation sigma_g > 1.");

  m.def("gamma_distribution", &lxmie::gammaDistribution,
    py::arg("effective_radius"), py::arg("effective_variance"),
    py::arg("nb_points") = lxmie::default_nb_distribution_points,
    py::arg("r_min") = 0.0, py::arg("r_max") = 0.0,
    "Gamma size distribution (Hansen & Travis 1974) with an effective radius and effective variance.");

  m.def("power_law_distribution", &lxmie::powerLawDistribution,
    py::arg("r_min"), py::arg("r_max"), py::arg("exponent"),
    py::arg("nb_points") = lxmie::default_nb_distribution_points,
    "Power-law size distribution n(r) ~ r^(-exponent) between r_min and r_max.");

  m.def("tabulated_distribution", &lxmie::tabulatedDistribution,
    py::arg("radius"), py::arg("dn_dr"),
    "Size distribution from tabulated values of n(r) = dN/dr at strictly increasing radii.");


  // -------------------------------------------------------------------------
  // MieDistributionResult
  // -------------------------------------------------------------------------

  py::class_<lxmie::MieDistributionResult>(m, "MieDistributionResult")
    .def_readonly("c_ext",                   &lxmie::MieDistributionResult::c_ext)
    .def_readonly("c_sca",                   &lxmie::MieDistributionResult::c_sca)
    .def_readonly("c_abs",                   &lxmie::MieDistributionResult::c_abs)
    .def_readonly("q_ext",                   &lxmie::MieDistributionResult::q_ext)
    .def_readonly("q_sca",                   &lxmie::MieDistributionResult::q_sca)
    .def_readonly("q_abs",                   &lxmie::MieDistributionResult::q_abs)
    .def_readonly("asymmetry_parameter",     &lxmie::MieDistributionResult::asymmetry_parameter)
    .def_readonly("geometric_cross_section", &lxmie::MieDistributionResult::geometric_cross_section)
    .def_readonly("effective_radius",        &lxmie::MieDistributionResult::effective_radius)
    .def("__repr__", [](const lxmie::MieDistributionResult& r) {
      return "MieDistributionResult(c_ext=" + std::to_string(r.c_ext)
           + ", c_sca=" + std::to_string(r.c_sca)
           + ", c_abs=" + std::to_string(r.c_abs)
           + ", asymmetry_parameter=" + std::to_string(r.asymmetry_parameter)
           + ", effective_radius=" + std::to_string(r.effective_radius) + ")";
    });


  // -------------------------------------------------------------------------
  // mie_distribution — optical properties integrated over a size distribution
  //   Overload 1: scalar n/wavelength → MieDistributionResult
  //   Overload 2: lists  n/wavelength → list[MieDistributionResult]
  // The wavelength has to be given in the same units as the particle radii.
  // -------------------------------------------------------------------------

  m.def("mie_distribution",
    [](double n_real, double n_imag, double wavelength, const lxmie::SizeDistribution& distribution) {
      py::gil_scoped_release release;
      return lxmie::MieDistribution(refractiveIndex(n_real, n_imag), wavelength, distribution);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"),
    "Compute cross sections and asymmetry parameter integrated over a size distribution (OpenMP-parallel over the particle sizes).");

  m.def("mie_distribution",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& wavelengths,
       const lxmie::SizeDistribution& distribution) {
      checkListSizes(n_real, n_imag, wavelengths);
      std::vector<lxmie::MieDistributionResult> results(n_real.size());
      parallelLoop(n_real.size(), [&](const std::size_t i) {
        results[i] = lxmie::MieDistribution(
          refractiveIndex(n_real[i], n_imag[i]), wavelengths[i], distribution);
      });
      return results;
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"),
    "Compute cross sections and asymmetry parameter integrated over a size distribution for a list of wavelengths (OpenMP-parallel).");


  // -------------------------------------------------------------------------
  // mie_distribution_legendre — additionally returns the Legendre series of
  // the size-averaged phase function
  //   Overload 1: scalar → (MieDistributionResult, list[float])
  //   Overload 2: lists  → (list[MieDistributionResult], list[list[float]])
  // -------------------------------------------------------------------------

  m.def("mie_distribution_legendre",
    [](double n_real, double n_imag, double wavelength,
       const lxmie::SizeDistribution& distribution, std::size_t nb_moments) {
      std::vector<double> moments;
      lxmie::MieDistributionResult r{};
      {
        py::gil_scoped_release release;
        r = lxmie::MieDistribution(
          refractiveIndex(n_real, n_imag), wavelength, distribution, nb_moments, moments);
      }
      return py::make_tuple(r, moments);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"), py::arg("nb_moments"),
    "Compute size-integrated cross sections plus the Legendre series of the size-averaged phase function.");

  m.def("mie_distribution_legendre",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& wavelengths,
       const lxmie::SizeDistribution& distribution,
       std::size_t nb_moments) {
      checkListSizes(n_real, n_imag, wavelengths);
      std::vector<lxmie::MieDistributionResult> results(n_real.size());
      std::vector<std::vector<double>> all_moments(n_real.size());
      parallelLoop(n_real.size(), [&](const std::size_t i) {
        results[i] = lxmie::MieDistribution(
          refractiveIndex(n_real[i], n_imag[i]), wavelengths[i], distribution, nb_moments, all_moments[i]);
      });
      return py::make_tuple(results, all_moments);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"), py::arg("nb_moments"),
    "Compute size-integrated cross sections plus the Legendre series of the size-averaged phase function for a list of wavelengths (OpenMP-parallel).");


  // -------------------------------------------------------------------------
  // mie_distribution_phase_function — additionally returns the size-averaged
  // phase function at the given angles (cosines of the scattering angles)
  //   Overload 1: scalar → (MieDistributionResult, list[float])
  //   Overload 2: lists  → (list[MieDistributionResult], list[list[float]])
  // -------------------------------------------------------------------------

  m.def("mie_distribution_phase_function",
    [](double n_real, double n_imag, double wavelength,
       const lxmie::SizeDistribution& distribution, const std::vector<double>& angles) {
      std::vector<double> phase_function;
      lxmie::MieDistributionResult r{};
      {
        py::gil_scoped_release release;
        r = lxmie::MieDistribution(
          refractiveIndex(n_real, n_imag), wavelength, distribution, angles, phase_function);
      }
      return py::make_tuple(r, phase_function);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"), py::arg("angles"),
    "Compute size-integrated cross sections plus the size-averaged phase function at the given angles.");

  m.def("mie_distribution_phase_function",
    [](const std::vector<double>& n_real,
       const std::vector<double>& n_imag,
       const std::vector<double>& wavelengths,
       const lxmie::SizeDistribution& distribution,
       const std::vector<double>& angles) {
      checkListSizes(n_real, n_imag, wavelengths);
      std::vector<lxmie::MieDistributionResult> results(n_real.size());
      std::vector<std::vector<double>> all_phase_functions(n_real.size());
      parallelLoop(n_real.size(), [&](const std::size_t i) {
        results[i] = lxmie::MieDistribution(
          refractiveIndex(n_real[i], n_imag[i]), wavelengths[i], distribution, angles, all_phase_functions[i]);
      });
      return py::make_tuple(results, all_phase_functions);
    },
    py::arg("n_real"), py::arg("n_imag"), py::arg("wavelength"), py::arg("distribution"), py::arg("angles"),
    "Compute size-integrated cross sections plus the size-averaged phase function at the given angles for a list of wavelengths (OpenMP-parallel).");
}

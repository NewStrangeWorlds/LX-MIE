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
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../mie/mie.h"

namespace py = pybind11;
using cd = std::complex<double>;


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
      return lxmie::Mie(cd(n_real, n_imag), x);
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
            results[i] = lxmie::Mie(cd(n_real[i], n_imag[i]), size_params[i]);
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
      auto r = lxmie::Mie(cd(n_real, n_imag), x, nb_moments, moments);
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
              cd(n_real[i], n_imag[i]), size_params[i], nb_moments, all_moments[i]);
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
      auto r = lxmie::Mie(cd(n_real, n_imag), x, angles, s1, s2);
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
              cd(n_real[i], n_imag[i]), size_params[i], angles, all_s1[i], all_s2[i]);
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
}

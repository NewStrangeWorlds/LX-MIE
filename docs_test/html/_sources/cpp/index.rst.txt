Using the library
=================

The Mie routines of LX-MIE are contained in a small C++17 library without
external dependencies. It consists of two header files and two source files
in the directory ``src/mie/``:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Files
     - Content
   * - ``mie.h``, ``mie.cpp``
     - Mie calculations for single particles
   * - ``size_distribution.h``, ``size_distribution.cpp``
     - Particle size distributions and size-integrated calculations

All functions and types are located in the namespace ``lxmie``. The complete
interface is documented in :doc:`api`.

Adding LX-MIE to a CMake project
--------------------------------

The directory ``src/mie/`` contains a ``CMakeLists.txt`` that defines the
static library target ``lxmie_lib``. It can be added to another CMake project
with ``add_subdirectory``:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.15)
   project(my_project LANGUAGES CXX)

   set(CMAKE_CXX_STANDARD 17)
   set(CMAKE_CXX_STANDARD_REQUIRED ON)

   # path to the LX-MIE repository
   add_subdirectory(lx-mie/src/mie lxmie_build)

   add_executable(my_code main.cpp)
   target_link_libraries(my_code PRIVATE lxmie_lib)

The target provides the include directory and, if available, links to OpenMP.
The headers are then included with

.. code-block:: cpp

   #include "mie.h"
   #include "size_distribution.h"

Alternatively, the two source files can simply be compiled together with the
other source files of a project, e.g.

.. code-block:: bash

   g++ -std=c++17 -O2 -fopenmp -I lx-mie/src/mie \
     main.cpp lx-mie/src/mie/mie.cpp lx-mie/src/mie/size_distribution.cpp -o my_code

A single particle
-----------------

.. code-block:: cpp

   #include <complex>
   #include <iostream>

   #include "mie.h"

   int main()
   {
     const double pi = 3.14159265358979323846;

     const double wavelength = 0.55;   //micron
     const double radius = 1.0;        //micron
     const double size_parameter = 2.0 * pi * radius / wavelength;

     const std::complex<double> refractive_index(1.5, 0.01);

     const lxmie::MieResult result = lxmie::Mie(refractive_index, size_parameter);

     std::cout << "Q_ext = " << result.q_ext << "\n"
               << "Q_sca = " << result.q_sca << "\n"
               << "Q_abs = " << result.q_abs << "\n"
               << "g     = " << result.asymmetry_parameter << "\n";

     return 0;
   }

Phase function and Legendre series
----------------------------------

The overloads of ``lxmie::Mie`` for the scattering amplitudes and the Legendre
series write their results to vectors that are passed by reference.

.. important::

   These output vectors have to be **allocated by the caller**. ``s1`` and
   ``s2`` need at least as many elements as the angular grid,
   ``legendre_moments`` at least ``nb_legendre_moments`` elements. The
   functions do not resize them.

.. code-block:: cpp

   #include <cmath>
   #include <complex>
   #include <vector>

   #include "mie.h"

   int main()
   {
     const std::complex<double> refractive_index(1.5, 0.01);
     const double size_parameter = 5.0;

     //cosines of the scattering angles: 0, 5, ..., 180 degrees
     std::vector<double> cos_angles;

     for (int i = 0; i <= 36; ++i)
       cos_angles.push_back(std::cos(i * 5.0 * 3.14159265358979323846 / 180.0));

     //scattering amplitudes and phase function
     std::vector<std::complex<double>> s1(cos_angles.size());
     std::vector<std::complex<double>> s2(cos_angles.size());

     const lxmie::MieResult result = lxmie::Mie(
       refractive_index, size_parameter, cos_angles, s1, s2);

     const std::vector<double> phase_function = lxmie::phaseFunction(
       s1, s2, size_parameter, result.q_sca);

     //Legendre series of the phase function
     const std::size_t nb_moments = 32;
     std::vector<double> legendre_moments(nb_moments);

     lxmie::Mie(refractive_index, size_parameter, nb_moments, legendre_moments);

     return 0;
   }

Size distributions
------------------

.. code-block:: cpp

   #include <complex>
   #include <iostream>
   #include <vector>

   #include "size_distribution.h"

   int main()
   {
     const std::complex<double> refractive_index(1.5, 0.01);
     const double wavelength = 0.55;   //micron

     //log-normal distribution, median radius 1 micron, sigma_g = 1.5
     const lxmie::SizeDistribution distribution = lxmie::logNormalDistribution(1.0, 1.5);

     const lxmie::MieDistributionResult result = lxmie::MieDistribution(
       refractive_index, wavelength, distribution);

     std::cout << "mean extinction cross section: " << result.c_ext << " micron^2\n"
               << "effective radius: " << result.effective_radius << " micron\n";

     //size-averaged phase function; the output vector is resized automatically
     const std::vector<double> cos_angles = {1.0, 0.5, 0.0, -0.5, -1.0};
     std::vector<double> phase_function;

     lxmie::MieDistribution(
       refractive_index, wavelength, distribution, cos_angles, phase_function);

     return 0;
   }

A size distribution can also be filled directly, e.g. with the size bins of a
cloud model:

.. code-block:: cpp

   lxmie::SizeDistribution distribution;

   distribution.radius = {0.01, 0.03, 0.1, 0.3, 1.0};      //bin radii
   distribution.weight = {2.0e3, 9.0e2, 2.5e2, 40.0, 3.0}; //number of particles per bin

Error handling
--------------

Errors are reported by exceptions:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Exception
     - Thrown for
   * - ``std::invalid_argument``
     - Invalid parameters of a size distribution or an invalid wavelength.
   * - ``std::runtime_error``
     - Failures within the Mie calculation, e.g. a continued fraction that
       did not converge. This is the case for a size parameter of zero.
   * - ``std::length_error``
     - A size parameter that is negative or not a number.

All messages start with ``LX-MIE:``.

Thread safety and parallelisation
---------------------------------

All functions of the library are thread-safe and can be called concurrently
from several threads, for example within an OpenMP-parallel loop over
wavelengths:

.. code-block:: cpp

   std::vector<lxmie::MieResult> results(nb_wavelengths);

   #pragma omp parallel for schedule(dynamic)
   for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(nb_wavelengths); ++i)
     results[i] = lxmie::Mie(refractive_index[i], size_parameter[i]);

.. note::

   Exceptions must not leave an OpenMP parallel region. If errors can occur
   within such a loop, they have to be caught inside the loop body.

The functions for single particles are not parallelised internally. If
compiled with OpenMP, ``lxmie::MieDistribution`` runs the loop over the
particle sizes in parallel, unless it is itself called from within a parallel
region. In that case, the loop over the sizes is executed serially, such that
the parallelisation of the calling code is not affected.

Internally, the routines keep their work arrays in thread-local storage and
reuse them in subsequent calls. This avoids repeated memory allocations, but
also means that the memory of the largest calculation of a thread remains
allocated until the thread ends (see :doc:`../background/performance`).

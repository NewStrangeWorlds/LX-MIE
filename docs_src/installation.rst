Installation
============

Python package
--------------

From PyPI
~~~~~~~~~

The easiest way to install the Python package is from
`PyPI <https://pypi.org/project/lxmie/>`_:

.. code-block:: bash

   pip install lxmie

Pre-built wheels are provided for CPython 3.9 – 3.13 on Linux (x86_64),
macOS (Apple Silicon), and Windows (x86_64). They include OpenMP support.
The package has no runtime dependencies.

.. note::

   Particle size distributions are available from version 1.2 onwards. The
   installed version can be checked with ``pip show lxmie``.

From the source code
~~~~~~~~~~~~~~~~~~~~

The package can also be compiled and installed from the repository:

.. code-block:: bash

   git clone https://github.com/newstrangeworlds/lx-mie.git
   cd lx-mie
   pip install .

This requires

- a C++17-compliant compiler (GCC ≥ 7, Clang ≥ 5, MSVC ≥ 2017),
- CMake ≥ 3.15,
- OpenMP (optional, but recommended).

The build system (``scikit-build-core`` and ``pybind11``) is downloaded
automatically by ``pip``.

To verify the installation, run

.. code-block:: bash

   python -c "import lxmie; print(lxmie.mie(1.5, 0.01, 5.0))"

which should print

.. code-block:: text

   MieResult(q_ext=3.818319, q_sca=3.554355, q_abs=0.263964, asymmetry_parameter=0.731372)


Command-line tool
-----------------

The command-line tool is built with CMake from the repository root
directory:

.. code-block:: bash

   git clone https://github.com/newstrangeworlds/lx-mie.git
   cd lx-mie
   cmake -B build
   cmake --build build

The executable ``lx_mie`` is placed in the root directory of the repository.
The requirements are the same as for the Python package: a C++17 compiler,
CMake ≥ 3.15, and optionally OpenMP.

By default, LX-MIE is compiled in ``Release`` mode. A different build type
can be chosen with the usual CMake option, e.g.
``cmake -B build -DCMAKE_BUILD_TYPE=Debug``.

To test the executable, run one of the :doc:`example parameter files
<cli/index>`:

.. code-block:: bash

   mkdir -p output
   ./lx_mie input/parameters.input


OpenMP
------

OpenMP is detected automatically by CMake. If it is available, the loops over
wavelengths, particles, or particle sizes are run in parallel (see
:doc:`background/performance`). Without OpenMP, LX-MIE compiles and runs
correctly in single-threaded mode.

On macOS, Apple Clang does not ship with OpenMP. It can be installed with
Homebrew (``brew install libomp``) and has to be made known to CMake, e.g.

.. code-block:: bash

   cmake -B build \
     -DOpenMP_CXX_FLAGS="-Xpreprocessor;-fopenmp;-I$(brew --prefix libomp)/include" \
     -DOpenMP_CXX_LIB_NAMES=omp \
     -DOpenMP_omp_LIBRARY=$(brew --prefix libomp)/lib/libomp.dylib


CMake options
-------------

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Option
     - Default
     - Description
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - Build type.
   * - ``BUILD_PYTHON_BINDINGS``
     - ``OFF``
     - Additionally builds the Python module with CMake. This requires an
       installation of `pybind11 <https://pybind11.readthedocs.io>`_ that can
       be found by CMake. The option is not needed for ``pip install .``.

The CMake project defines three targets: the static library ``lxmie_lib``
with the Mie routines, the executable ``lx_mie``, and, optionally, the Python
module ``lxmie``. Using the library in other C++ projects is described in
:doc:`cpp/index`.

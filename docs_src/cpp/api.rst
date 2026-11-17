API reference
=============

.. default-domain:: cpp

All functions and types are declared in the namespace ``lxmie``. The
definitions of the physical quantities are given in :doc:`../conventions`.

.. cpp:namespace:: lxmie

Single particles
----------------

Header: ``mie.h``

.. cpp:struct:: MieResult

   Result of a Mie calculation for a single particle.

   .. cpp:member:: double q_ext

      Extinction efficiency.

   .. cpp:member:: double q_sca

      Scattering efficiency.

   .. cpp:member:: double q_abs

      Absorption efficiency, calculated as the difference of the extinction
      and scattering efficiencies.

   .. cpp:member:: double asymmetry_parameter

      Asymmetry parameter :math:`g`.


.. cpp:function:: MieResult Mie(std::complex<double> refractive_index, double size_parameter)

   Calculates the Mie efficiencies and the asymmetry parameter.

   :param refractive_index: Complex refractive index of the particle. The
      sign of the imaginary part does not matter.
   :param size_parameter: Size parameter :math:`x = 2 \pi r / \lambda`.
   :throws std\:\:runtime_error: if the calculation fails.


.. cpp:function:: MieResult Mie(std::complex<double> refractive_index, double size_parameter, std::size_t nb_legendre_moments, std::vector<double>& legendre_moments)

   Additionally calculates the coefficients of the
   :ref:`Legendre series <legendre-series>` of the phase function.

   :param refractive_index: Complex refractive index of the particle.
   :param size_parameter: Size parameter.
   :param nb_legendre_moments: Number of Legendre series coefficients.
   :param legendre_moments: Output: the coefficients
      :math:`\chi_0, \chi_1, \ldots`. The vector has to contain at least
      ``nb_legendre_moments`` elements. It is not resized by the function.
   :throws std\:\:runtime_error: if the calculation fails.


.. cpp:function:: MieResult Mie(std::complex<double> refractive_index, double size_parameter, const std::vector<double>& angular_grid, std::vector<std::complex<double>>& s1, std::vector<std::complex<double>>& s2)

   Additionally calculates the scattering amplitudes :math:`S_1` and
   :math:`S_2` on a grid of scattering angles.

   :param refractive_index: Complex refractive index of the particle.
   :param size_parameter: Size parameter.
   :param angular_grid: **Cosines** of the scattering angles.
   :param s1: Output: scattering amplitude :math:`S_1` at each angle.
   :param s2: Output: scattering amplitude :math:`S_2` at each angle.
   :throws std\:\:runtime_error: if the calculation fails.

   The vectors ``s1`` and ``s2`` have to contain at least as many elements as
   ``angular_grid``. They are not resized by the function.


.. cpp:function:: double phaseFunction(std::complex<double> s1, std::complex<double> s2, double size_parameter, double q_sca)

   Returns the :ref:`normalised phase function <phase-function-normalisation>`
   for a single scattering angle,

   .. math::

      p = \frac{2}{x^2 \, Q_\mathrm{sca}} \left( |S_1|^2 + |S_2|^2 \right) \ .

   :param s1: Scattering amplitude :math:`S_1`.
   :param s2: Scattering amplitude :math:`S_2`.
   :param size_parameter: Size parameter of the particle.
   :param q_sca: Scattering efficiency of the particle.


.. cpp:function:: std::vector<double> phaseFunction(const std::vector<std::complex<double>>& s1, const std::vector<std::complex<double>>& s2, double size_parameter, double q_sca)

   Returns the normalised phase function for a set of scattering angles.

   :param s1: Scattering amplitudes :math:`S_1`.
   :param s2: Scattering amplitudes :math:`S_2`, with at least as many
      elements as ``s1``.
   :param size_parameter: Size parameter of the particle.
   :param q_sca: Scattering efficiency of the particle.


Size distributions
------------------

Header: ``size_distribution.h``

.. cpp:var:: constexpr std::size_t default_nb_distribution_points = 500

   Default number of radius grid points of the analytical size distributions.


.. cpp:struct:: SizeDistribution

   A particle size distribution, discretised into a set of quadrature nodes.
   The struct can be filled directly, e.g. with the size bins of a cloud
   model, or created with one of the functions below.

   .. cpp:member:: std::vector<double> radius

      Particle radii.

   .. cpp:member:: std::vector<double> weight

      Number of particles represented by each node, i.e.
      :math:`n(r_i) \, \Delta r_i` for a continuous distribution
      :math:`n(r) = \mathrm{d}N/\mathrm{d}r`. The weights do not need to be
      normalised. All size-integrated results are normalised to the total
      number of particles, i.e. the sum of all weights.


.. cpp:function:: SizeDistribution logNormalDistribution(double median_radius, double sigma_g, std::size_t nb_points = default_nb_distribution_points, double r_min = 0.0, double r_max = 0.0)

   Log-normal distribution with the median radius :math:`r_m` and the
   geometric standard deviation :math:`\sigma_g`,

   .. math::

      n(r) \propto \frac{1}{r} \exp \left( - \frac{\ln^2 (r / r_m)}{2 \ln^2 \sigma_g} \right) \ .

   :param median_radius: Median radius, has to be positive.
   :param sigma_g: Geometric standard deviation, has to be larger than 1.
   :param nb_points: Number of radius grid points (at least 2).
   :param r_min: Smallest radius. If zero, the integration limit is chosen
      automatically.
   :param r_max: Largest radius. If zero, the integration limit is chosen
      automatically.
   :throws std\:\:invalid_argument: for invalid parameters.


.. cpp:function:: SizeDistribution gammaDistribution(double effective_radius, double effective_variance, std::size_t nb_points = default_nb_distribution_points, double r_min = 0.0, double r_max = 0.0)

   Gamma distribution following Hansen & Travis (1974) with the effective
   radius :math:`a` and the effective variance :math:`b`,

   .. math::

      n(r) \propto r^{(1-3b)/b} \exp \left( - \frac{r}{ab} \right) \ .

   :param effective_radius: Effective radius, has to be positive.
   :param effective_variance: Effective variance, between 0 and 0.5.
   :param nb_points: Number of radius grid points (at least 2).
   :param r_min: Smallest radius. If zero, the integration limit is chosen
      automatically.
   :param r_max: Largest radius. If zero, the integration limit is chosen
      automatically.
   :throws std\:\:invalid_argument: for invalid parameters.


.. cpp:function:: SizeDistribution powerLawDistribution(double r_min, double r_max, double exponent, std::size_t nb_points = default_nb_distribution_points)

   Power-law distribution :math:`n(r) \propto r^{-p}` between ``r_min`` and
   ``r_max``.

   :param r_min: Smallest radius, has to be positive.
   :param r_max: Largest radius, has to be larger than ``r_min``.
   :param exponent: Exponent :math:`p`.
   :param nb_points: Number of radius grid points (at least 2).
   :throws std\:\:invalid_argument: for invalid parameters.


.. cpp:function:: SizeDistribution tabulatedDistribution(const std::vector<double>& radius, const std::vector<double>& dn_dr)

   Tabulated distribution :math:`n(r) = \mathrm{d}N/\mathrm{d}r`. The weights
   are obtained from the trapezoidal rule.

   :param radius: Particle radii, strictly increasing (at least two).
   :param dn_dr: Values of the distribution at these radii.
   :throws std\:\:invalid_argument: for vectors of different sizes, fewer
      than two radii, or radii that are not strictly increasing.


Size-integrated calculations
----------------------------

Header: ``size_distribution.h``

.. cpp:struct:: MieDistributionResult

   Size-integrated optical properties, see :ref:`size-integrated-quantities`.
   Cross sections are mean values per particle, in units of the squared
   length unit that is used for the wavelength and the particle radii.

   .. cpp:member:: double c_ext

      Mean extinction cross section.

   .. cpp:member:: double c_sca

      Mean scattering cross section.

   .. cpp:member:: double c_abs

      Mean absorption cross section.

   .. cpp:member:: double q_ext

      Effective extinction efficiency, i.e. the mean extinction cross section
      divided by the mean geometric cross section.

   .. cpp:member:: double q_sca

      Effective scattering efficiency.

   .. cpp:member:: double q_abs

      Effective absorption efficiency.

   .. cpp:member:: double asymmetry_parameter

      Asymmetry parameter, weighted by the scattering cross sections.

   .. cpp:member:: double geometric_cross_section

      Mean geometric cross section :math:`\langle \pi r^2 \rangle`.

   .. cpp:member:: double effective_radius

      Effective radius :math:`\langle r^3 \rangle / \langle r^2 \rangle`.


.. cpp:function:: MieDistributionResult MieDistribution(std::complex<double> refractive_index, double wavelength, const SizeDistribution& distribution)

   Mie calculation integrated over a particle size distribution.

   :param refractive_index: Complex refractive index of the particles. The
      sign of the imaginary part does not matter.
   :param wavelength: Wavelength, in the same units as the particle radii.
   :param distribution: The size distribution.
   :throws std\:\:invalid_argument: if the wavelength is not positive or the
      distribution is invalid (different numbers of radii and weights, no
      radii, radii that are not positive, negative weights, or a sum of the
      weights that is not positive).
   :throws std\:\:runtime_error: if one of the Mie calculations fails.

   Nodes with a weight of zero are skipped.


.. cpp:function:: MieDistributionResult MieDistribution(std::complex<double> refractive_index, double wavelength, const SizeDistribution& distribution, std::size_t nb_legendre_moments, std::vector<double>& legendre_moments)

   Additionally returns the Legendre series of the size-averaged phase
   function.

   :param refractive_index: Complex refractive index of the particles.
   :param wavelength: Wavelength, in the same units as the particle radii.
   :param distribution: The size distribution.
   :param nb_legendre_moments: Number of Legendre series coefficients.
   :param legendre_moments: Output: the coefficients of the Legendre series.
      The vector is resized to ``nb_legendre_moments`` elements.


.. cpp:function:: MieDistributionResult MieDistribution(std::complex<double> refractive_index, double wavelength, const SizeDistribution& distribution, const std::vector<double>& angular_grid, std::vector<double>& phase_function)

   Additionally returns the size-averaged phase function.

   :param refractive_index: Complex refractive index of the particles.
   :param wavelength: Wavelength, in the same units as the particle radii.
   :param distribution: The size distribution.
   :param angular_grid: **Cosines** of the scattering angles.
   :param phase_function: Output: the size-averaged phase function at each
      angle. The vector is resized to the size of ``angular_grid``.

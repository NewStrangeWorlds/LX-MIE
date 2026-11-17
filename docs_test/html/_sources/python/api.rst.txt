API reference
=============

.. py:module:: lxmie

This page documents all functions and classes of the ``lxmie`` module. The
general rules for scalar and batch calls are described in :doc:`index`, the
definitions of the physical quantities in :doc:`../conventions`.

In the signatures below, ``float | list[float]`` denotes an argument that is
either a scalar or a list. If lists are used, all arguments of this kind have
to be lists of the same length and the function returns lists of results.


Single particles
----------------

.. py:function:: mie(n_real, n_imag, size_parameter)

   Calculates the Mie efficiencies and the asymmetry parameter.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param size_parameter: Size parameter :math:`x = 2 \pi r / \lambda`.
   :type size_parameter: float | list[float]
   :returns: The efficiencies and the asymmetry parameter.
   :rtype: MieResult | list[MieResult]

   .. code-block:: python

      result = lxmie.mie(1.5, 0.01, 5.0)
      results = lxmie.mie([1.5, 1.6], [0.01, 0.02], [5.0, 10.0])


.. py:function:: mie_legendre(n_real, n_imag, size_parameter, nb_moments)

   Calculates the Mie efficiencies, the asymmetry parameter, and the first
   ``nb_moments`` coefficients of the
   :ref:`Legendre series <legendre-series>` of the phase function.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param size_parameter: Size parameter :math:`x = 2 \pi r / \lambda`.
   :type size_parameter: float | list[float]
   :param int nb_moments: Number of Legendre series coefficients, starting
      with :math:`\chi_0`.
   :returns: A tuple ``(result, moments)``. ``moments`` is the list of
      Legendre coefficients :math:`\chi_0 \ldots \chi_{\mathrm{nb\_moments}-1}`.
      For list arguments, both elements of the tuple are lists with one entry
      per particle.
   :rtype: tuple[MieResult, list[float]] | tuple[list[MieResult], list[list[float]]]

   The first coefficient is always 1, the second one is equal to the
   asymmetry parameter.

   .. code-block:: python

      result, moments = lxmie.mie_legendre(1.5, 0.01, 5.0, 16)

      results, all_moments = lxmie.mie_legendre([1.5, 1.6], [0.01, 0.02], [5.0, 10.0], 16)
      # all_moments[i][l]: coefficient l of particle i


.. py:function:: mie_amplitudes(n_real, n_imag, size_parameter, angles)

   Calculates the Mie efficiencies, the asymmetry parameter, and the complex
   scattering amplitudes :math:`S_1` and :math:`S_2` at the given angles.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param size_parameter: Size parameter :math:`x = 2 \pi r / \lambda`.
   :type size_parameter: float | list[float]
   :param list[float] angles: **Cosines** of the scattering angles. The same
      angles are used for all particles of a batch call.
   :returns: A tuple ``(result, s1, s2)``. ``s1`` and ``s2`` are lists with
      the scattering amplitudes at each angle. For list arguments, all three
      elements of the tuple are lists with one entry per particle.
   :rtype: tuple[MieResult, list[complex], list[complex]] | tuple[list[MieResult], list[list[complex]], list[list[complex]]]

   The scattering amplitudes can be converted to the phase function with
   :func:`phase_function`.

   .. code-block:: python

      import math

      cos_angles = [math.cos(math.radians(theta)) for theta in range(0, 181, 5)]

      result, s1, s2 = lxmie.mie_amplitudes(1.5, 0.01, 5.0, cos_angles)

   .. note::

      Due to the internal sign convention of the refractive index
      (:math:`m = n - \mathrm{i}k`), the scattering amplitudes are the
      complex conjugates of those obtained with the convention
      :math:`m = n + \mathrm{i}k` that is used by, for example, Bohren &
      Huffman (1983). Intensities, such as :math:`|S_1|^2` and
      :math:`|S_2|^2`, are not affected.


.. py:function:: phase_function(s1, s2, size_parameter, q_sca)

   Calculates the :ref:`normalised phase function
   <phase-function-normalisation>` for unpolarised incident light from the
   scattering amplitudes,

   .. math::

      p = \frac{2}{x^2 \, Q_\mathrm{sca}} \left( |S_1|^2 + |S_2|^2 \right) \ .

   :param s1: Scattering amplitude(s) :math:`S_1`.
   :type s1: complex | list[complex]
   :param s2: Scattering amplitude(s) :math:`S_2`.
   :type s2: complex | list[complex]
   :param float size_parameter: Size parameter of the particle.
   :param float q_sca: Scattering efficiency of the particle.
   :returns: The phase function at the angle(s) of the scattering amplitudes.
   :rtype: float | list[float]

   In contrast to the other functions, the list version of this function
   refers to a list of angles for a *single* particle. For the result of a
   batch call of :func:`mie_amplitudes`, the function has to be called for
   each particle:

   .. code-block:: python

      results, all_s1, all_s2 = lxmie.mie_amplitudes(n_real, n_imag, size_parameters, cos_angles)

      phase_functions = [
          lxmie.phase_function(s1, s2, x, r.q_sca)
          for s1, s2, x, r in zip(all_s1, all_s2, size_parameters, results)
      ]


.. py:class:: MieResult

   Result of a Mie calculation for a single particle. All attributes are
   read-only.

   .. py:attribute:: q_ext
      :type: float

      Extinction efficiency.

   .. py:attribute:: q_sca
      :type: float

      Scattering efficiency.

   .. py:attribute:: q_abs
      :type: float

      Absorption efficiency, :math:`Q_\mathrm{ext} - Q_\mathrm{sca}`.

   .. py:attribute:: asymmetry_parameter
      :type: float

      Asymmetry parameter :math:`g`.

   .. note::

      The printed representation of the result objects, e.g.
      ``MieResult(q_ext=3.818319, ...)``, is rounded to six decimal places.
      The attributes contain the values with full double precision.


Size distributions
------------------

A size distribution is represented by a :class:`SizeDistribution` object. It
can be created with one of the following functions or directly from lists of
radii and weights. The definitions of the distributions and the choice of the
radius grids are discussed in :doc:`../background/size_distributions`.

.. py:function:: lognormal_distribution(median_radius, sigma_g, nb_points=500, r_min=0, r_max=0)

   Log-normal distribution,

   .. math::

      n(r) \propto \frac{1}{r} \exp \left( - \frac{\ln^2 (r / r_m)}{2 \ln^2 \sigma_g} \right) \ .

   :param float median_radius: Median radius :math:`r_m > 0`.
   :param float sigma_g: Geometric standard deviation :math:`\sigma_g > 1`.
   :param int nb_points: Number of radii used for the integration over the
      distribution (at least 2).
   :param float r_min: Smallest radius of the grid. If zero, the limit is
      chosen automatically.
   :param float r_max: Largest radius of the grid. If zero, the limit is
      chosen automatically.
   :rtype: SizeDistribution
   :raises ValueError: for invalid parameters.


.. py:function:: gamma_distribution(effective_radius, effective_variance, nb_points=500, r_min=0, r_max=0)

   Gamma distribution of Hansen & Travis (1974),

   .. math::

      n(r) \propto r^{(1-3b)/b} \exp \left( - \frac{r}{ab} \right) \ ,

   with the effective radius :math:`a` and the effective variance :math:`b`.

   :param float effective_radius: Effective radius :math:`a > 0`.
   :param float effective_variance: Effective variance :math:`0 < b < 0.5`.
   :param int nb_points: Number of radii used for the integration over the
      distribution (at least 2).
   :param float r_min: Smallest radius of the grid. If zero, the limit is
      chosen automatically.
   :param float r_max: Largest radius of the grid. If zero, the limit is
      chosen automatically.
   :rtype: SizeDistribution
   :raises ValueError: for invalid parameters.


.. py:function:: power_law_distribution(r_min, r_max, exponent, nb_points=500)

   Power-law distribution :math:`n(r) \propto r^{-p}` between ``r_min`` and
   ``r_max``.

   :param float r_min: Smallest radius, :math:`r_\mathrm{min} > 0`.
   :param float r_max: Largest radius,
      :math:`r_\mathrm{max} > r_\mathrm{min}`.
   :param float exponent: Exponent :math:`p`. Note the sign: a positive
      value describes a distribution that decreases with the radius.
   :param int nb_points: Number of radii used for the integration over the
      distribution (at least 2).
   :rtype: SizeDistribution
   :raises ValueError: for invalid parameters.


.. py:function:: tabulated_distribution(radius, dn_dr)

   Distribution from tabulated values of :math:`n(r) = \mathrm{d}N/\mathrm{d}r`.
   The integration weights are calculated with the trapezoidal rule.

   :param list[float] radius: Particle radii, strictly increasing (at least
      two values).
   :param list[float] dn_dr: Values of the distribution at these radii, in
      arbitrary units.
   :rtype: SizeDistribution
   :raises ValueError: for lists of different lengths, fewer than two radii,
      or radii that are not strictly increasing.


.. py:class:: SizeDistribution(radius, weight)

   A particle size distribution, discretised into a set of radii and the
   number of particles represented by each radius.

   The constructor directly uses the given radii and weights. This is the
   appropriate choice for discrete size bins, e.g. from a cloud model, with
   the number (density) of particles in each bin as weights. In contrast to
   :func:`tabulated_distribution`, the weights are *not* values of
   :math:`\mathrm{d}N/\mathrm{d}r`, but particle numbers.

   :param list[float] radius: Particle radii. They have to be positive, but
      do not need to be sorted.
   :param list[float] weight: Number of particles for each radius. The
      weights have to be non-negative and do not need to be normalised.
   :raises ValueError: if the two lists have different lengths.

   .. py:attribute:: radius
      :type: list[float]

      The radii of the distribution.

   .. py:attribute:: weight
      :type: list[float]

      The weights of the distribution. For a distribution created from a
      continuous function :math:`n(r)`, these are the quadrature weights
      :math:`n(r_i) \, \Delta r_i`.

   ``len(distribution)`` returns the number of radii.

   .. note::

      ``radius`` and ``weight`` return copies of the internal data. To change
      a distribution, a complete new list has to be assigned, e.g.
      ``distribution.weight = new_weights``. Changing single elements, such
      as ``distribution.weight[0] = 0``, has no effect.

   The radii and weights are validated when the distribution is used in a
   calculation, see :func:`mie_distribution`.


Size-integrated calculations
----------------------------

The following functions integrate the optical properties over a size
distribution. The calculated quantities are defined in
:ref:`size-integrated-quantities`.

.. py:function:: mie_distribution(n_real, n_imag, wavelength, distribution)

   Calculates the cross sections, efficiencies, and the asymmetry parameter
   integrated over a size distribution.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param wavelength: Wavelength, in the same units as the radii of the
      distribution.
   :type wavelength: float | list[float]
   :param SizeDistribution distribution: The size distribution. The same
      distribution is used for all wavelengths of a batch call.
   :returns: The size-integrated optical properties.
   :rtype: MieDistributionResult | list[MieDistributionResult]
   :raises ValueError: if the lists have different lengths, the wavelength
      is not positive, or the distribution is invalid (different numbers of
      radii and weights, no radii, radii that are not positive, negative
      weights, or a sum of the weights that is not positive).

   .. code-block:: python

      distribution = lxmie.lognormal_distribution(1.0, 1.5)   # radii in micron

      result = lxmie.mie_distribution(1.5, 0.01, 0.55, distribution)
      results = lxmie.mie_distribution([1.5, 1.6], [0.01, 0.02], [0.55, 1.0], distribution)


.. py:function:: mie_distribution_legendre(n_real, n_imag, wavelength, distribution, nb_moments)

   Like :func:`mie_distribution`, but additionally returns the first
   ``nb_moments`` coefficients of the Legendre series of the size-averaged
   phase function. The coefficients have the same
   :ref:`normalisation <legendre-series>` as those of :func:`mie_legendre`.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param wavelength: Wavelength, in the same units as the radii of the
      distribution.
   :type wavelength: float | list[float]
   :param SizeDistribution distribution: The size distribution.
   :param int nb_moments: Number of Legendre series coefficients.
   :returns: A tuple ``(result, moments)``. For list arguments, both elements
      are lists with one entry per wavelength.
   :rtype: tuple[MieDistributionResult, list[float]] | tuple[list[MieDistributionResult], list[list[float]]]
   :raises ValueError: see :func:`mie_distribution`.


.. py:function:: mie_distribution_phase_function(n_real, n_imag, wavelength, distribution, angles)

   Like :func:`mie_distribution`, but additionally returns the size-averaged
   phase function at the given angles. The phase function has the same
   :ref:`normalisation <phase-function-normalisation>` as the one returned by
   :func:`phase_function`.

   :param n_real: Real part of the refractive index.
   :type n_real: float | list[float]
   :param n_imag: Imaginary part of the refractive index. The sign does not
      matter.
   :type n_imag: float | list[float]
   :param wavelength: Wavelength, in the same units as the radii of the
      distribution.
   :type wavelength: float | list[float]
   :param SizeDistribution distribution: The size distribution.
   :param list[float] angles: **Cosines** of the scattering angles.
   :returns: A tuple ``(result, phase_function)``. For list arguments, both
      elements are lists with one entry per wavelength.
   :rtype: tuple[MieDistributionResult, list[float]] | tuple[list[MieDistributionResult], list[list[float]]]
   :raises ValueError: see :func:`mie_distribution`.

   .. code-block:: python

      import math

      cos_angles = [math.cos(math.radians(theta)) for theta in range(0, 181, 15)]

      result, phase_function = lxmie.mie_distribution_phase_function(
          1.5, 0.01, 0.55, distribution, cos_angles)


.. py:class:: MieDistributionResult

   Result of a calculation for a size distribution. All attributes are
   read-only. The length unit is the one used for the wavelength and the radii
   of the distribution.

   .. py:attribute:: c_ext
      :type: float

      Mean extinction cross section per particle (length unit squared).

   .. py:attribute:: c_sca
      :type: float

      Mean scattering cross section per particle (length unit squared).

   .. py:attribute:: c_abs
      :type: float

      Mean absorption cross section per particle (length unit squared).

   .. py:attribute:: q_ext
      :type: float

      Effective extinction efficiency, i.e. the mean extinction cross
      section divided by the mean geometric cross section.

   .. py:attribute:: q_sca
      :type: float

      Effective scattering efficiency.

   .. py:attribute:: q_abs
      :type: float

      Effective absorption efficiency.

   .. py:attribute:: asymmetry_parameter
      :type: float

      Asymmetry parameter of the distribution, weighted by the scattering
      cross sections.

   .. py:attribute:: geometric_cross_section
      :type: float

      Mean geometric cross section :math:`\langle \pi r^2 \rangle` (length
      unit squared).

   .. py:attribute:: effective_radius
      :type: float

      Effective radius :math:`\langle r^3 \rangle / \langle r^2 \rangle`
      (length unit).

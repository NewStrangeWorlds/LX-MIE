Conventions and units
=====================

This page summarises the definitions that are used by all three interfaces of
LX-MIE. The underlying equations are given in :doc:`background/theory`.

Size parameter
--------------

The size parameter of a spherical particle with the radius :math:`r` at the
wavelength :math:`\lambda` is

.. math::

   x = \frac{2 \pi r}{\lambda} \ .

It is dimensionless. The Mie efficiencies only depend on the size parameter
and the refractive index. The single-particle functions of the Python package
and the C++ library therefore take the size parameter as their argument,
rather than the radius and the wavelength.

The wavelength refers to the medium surrounding the particle and the
refractive index is the one of the particle relative to that medium. For
particles in a gas, this is in practice the vacuum wavelength and the
refractive index of the particle material.

Refractive index
----------------

The complex refractive index is given by its real part :math:`n` and its
imaginary part :math:`k`.

**The sign of the imaginary part does not matter.** Internally, LX-MIE uses
the convention :math:`m = n - \mathrm{i} k` and always converts the imaginary
part to a negative value. Both sign conventions, :math:`n + \mathrm{i} k` and
:math:`n - \mathrm{i} k`, therefore give the same, physically correct result
for an absorbing particle. The values of :math:`k` from the files of the
:doc:`compilation <background/refractive_indices>` can be used directly.

In the Python package, the refractive index is always passed as two separate
real numbers (``n_real``, ``n_imag``). The C++ library uses
``std::complex<double>``.

Efficiencies and cross sections
-------------------------------

LX-MIE calculates the extinction, scattering, and absorption efficiencies
:math:`Q_\mathrm{ext}`, :math:`Q_\mathrm{sca}`, and :math:`Q_\mathrm{abs}`,
with

.. math::

   Q_\mathrm{ext} = Q_\mathrm{sca} + Q_\mathrm{abs} \ .

The corresponding cross sections are obtained by multiplying the efficiencies
with the geometric cross section of the particle:

.. math::

   C = Q \, \pi r^2 \ .

The absorption efficiency is calculated as the difference of the extinction
and scattering efficiencies. For non-absorbing particles (:math:`k = 0`), it
can therefore deviate from zero by the numerical round-off error, including
very small negative values.

Other frequently used quantities directly follow from the results, for
example the single scattering albedo
:math:`\omega_0 = Q_\mathrm{sca} / Q_\mathrm{ext}`.

Asymmetry parameter
-------------------

The asymmetry parameter :math:`g` is the mean cosine of the scattering angle,

.. math::

   g = \langle \cos\theta \rangle
     = \frac{1}{2} \int_{-1}^{1} p(\mu) \, \mu \, \mathrm{d}\mu \ ,

with :math:`\mu = \cos\theta`. It is :math:`g = 0` for isotropic or symmetric
scattering (e.g. Rayleigh scattering) and approaches :math:`g = 1` for pure
forward scattering.

Scattering angles
-----------------

Wherever angles are passed to or returned by LX-MIE, they are given as the
**cosines of the scattering angles**, :math:`\mu = \cos\theta`, and not as
angles in degrees or radians. Forward scattering corresponds to
:math:`\mu = 1`, backward scattering to :math:`\mu = -1`.

.. _phase-function-normalisation:

Phase function
--------------

The phase function :math:`p(\mu)` for unpolarised incident light is normalised
such that

.. math::

   \frac{1}{2} \int_{-1}^{1} p(\mu) \, \mathrm{d}\mu = 1
   \qquad \Longleftrightarrow \qquad
   \int_{4\pi} p \, \mathrm{d}\Omega = 4 \pi \ .

Isotropic scattering therefore corresponds to :math:`p(\mu) = 1`.

.. _legendre-series:

Legendre series
---------------

The Legendre series coefficients :math:`\chi_l` returned by LX-MIE are defined
by

.. math::

   p(\mu) = \sum_{l=0}^{\infty} \left( 2l + 1 \right) \chi_l \, P_l(\mu) \ ,
   \qquad
   \chi_l = \frac{1}{2} \int_{-1}^{1} p(\mu) \, P_l(\mu) \, \mathrm{d}\mu \ ,

where :math:`P_l` are the Legendre polynomials. With the normalisation of the
phase function given above, the first two coefficients are

.. math::

   \chi_0 = 1 \ , \qquad \chi_1 = g \ .

This is the convention used by, for example, the discrete ordinate radiative
transfer code DISORT. Other codes include the factor :math:`(2l+1)` in the
coefficients. In that case, the coefficients from LX-MIE have to be multiplied
by :math:`(2l+1)`.

The series of a particle with the size parameter :math:`x` contains about
:math:`2x` non-negligible terms. Coefficients beyond the end of the series
are returned as zero.

Units
-----

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Interface
     - Units
   * - Python, C++: single particles
     - Only dimensionless quantities are used (size parameter, refractive
       index, efficiencies).
   * - Python, C++: size distributions
     - Radii and wavelengths can be given in any length unit, as long as it
       is the same for both. Cross sections are returned in the square of
       that unit, e.g. in |micron|\ :sup:`2` for radii and wavelengths in
       |micron|.
   * - Command-line tool
     - Radii and wavelengths in |micron|, cross sections in
       cm\ :sup:`2`.

Limitations
-----------

LX-MIE treats homogeneous, isotropic, spherical particles. Coated spheres,
non-spherical particles, and magnetic materials are not supported. The
phase function is the one for unpolarised incident light; the complex
scattering amplitudes :math:`S_1` and :math:`S_2`, from which all elements of
the scattering matrix can be constructed, are available via
:func:`lxmie.mie_amplitudes`.

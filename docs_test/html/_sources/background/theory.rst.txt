Theory and numerical method
===========================

This page summarises the equations that are evaluated by LX-MIE and the
numerical methods that are used to do so. A detailed description is given
in Kitzmann & Heng (2018) [KH18]_. Comprehensive introductions to Mie theory
can be found in the textbooks by Bohren & Huffman [BH83]_ and van de Hulst
[vdH57]_.

Mie theory
----------

Mie theory [Mie08]_ describes the scattering and absorption of a plane
electromagnetic wave by a homogeneous sphere. The solution only depends on
two quantities: the complex refractive index :math:`m` of the sphere relative
to the surrounding medium and the size parameter

.. math::

   x = \frac{2 \pi r}{\lambda} \ ,

where :math:`r` is the radius of the sphere and :math:`\lambda` the wavelength
in the surrounding medium.

All optical properties are expressed as infinite series in terms of the Mie
coefficients :math:`a_n` and :math:`b_n`. The extinction and scattering
efficiencies are given by

.. math::

   Q_\mathrm{ext} = \frac{2}{x^2} \sum_{n=1}^{\infty} \left( 2n + 1 \right)
                    \mathrm{Re} \left( a_n + b_n \right) \ ,

.. math::

   Q_\mathrm{sca} = \frac{2}{x^2} \sum_{n=1}^{\infty} \left( 2n + 1 \right)
                    \left( |a_n|^2 + |b_n|^2 \right) \ .

The absorption efficiency follows from
:math:`Q_\mathrm{abs} = Q_\mathrm{ext} - Q_\mathrm{sca}`.

The asymmetry parameter is calculated from

.. math::

   g = \frac{4}{x^2 Q_\mathrm{sca}} \sum_{n=1}^{\infty} \left[
       \frac{n \left( n + 2 \right)}{n + 1}
       \mathrm{Re} \left( a_n a_{n+1}^* + b_n b_{n+1}^* \right)
       + \frac{2n + 1}{n \left( n + 1 \right)}
       \mathrm{Re} \left( a_n b_n^* \right) \right] \ ,

where the asterisk denotes the complex conjugate.

Scattering amplitudes and phase function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The angular distribution of the scattered light is described by the two
complex scattering amplitudes

.. math::

   S_1(\mu) = \sum_{n=1}^{\infty} \frac{2n + 1}{n \left( n + 1 \right)}
              \left[ a_n \pi_n(\mu) + b_n \tau_n(\mu) \right] \ ,

.. math::

   S_2(\mu) = \sum_{n=1}^{\infty} \frac{2n + 1}{n \left( n + 1 \right)}
              \left[ a_n \tau_n(\mu) + b_n \pi_n(\mu) \right] \ ,

with the cosine of the scattering angle :math:`\mu = \cos\theta` and the
angular functions :math:`\pi_n` and :math:`\tau_n`, which are related to the
Legendre polynomials and their derivatives. The phase function for
unpolarised incident light is

.. math::

   p(\mu) = \frac{2}{x^2 Q_\mathrm{sca}} \left( |S_1(\mu)|^2 + |S_2(\mu)|^2 \right) \ .

It is :ref:`normalised <phase-function-normalisation>` such that its average
over all directions is unity.

Numerical method
----------------

The direct evaluation of the Mie coefficients from their definition in terms
of spherical Bessel functions is numerically unstable, in particular for
large size parameters and strongly absorbing particles. LX-MIE therefore
expresses the Mie coefficients by ratios of Riccati-Bessel functions, which
can be calculated by stable recurrence relations.

Number of terms
~~~~~~~~~~~~~~~

The Mie series are truncated after

.. math::

   N = x + 4.3 \, x^{1/3} + 2

terms (rounded down to an integer). The computing time and the memory
required for a calculation therefore increase linearly with the size
parameter (see :doc:`performance`).

Mie coefficients
~~~~~~~~~~~~~~~~

The Mie coefficients are calculated from

.. math::

   a_n = C_n(x) \, \frac{A_n(mx) / m - A_n(x)}{A_n(mx) / m - B_n(x)} \ ,
   \qquad
   b_n = C_n(x) \, \frac{m \, A_n(mx) - A_n(x)}{m \, A_n(mx) - B_n(x)} \ ,

where :math:`A_n` and :math:`B_n` are the logarithmic derivatives of the
Riccati-Bessel functions :math:`\psi_n` and :math:`\zeta_n`, respectively,
and :math:`C_n` is their ratio:

.. math::

   A_n(z) = \frac{\psi_n'(z)}{\psi_n(z)} \ , \qquad
   B_n(x) = \frac{\zeta_n'(x)}{\zeta_n(x)} \ , \qquad
   C_n(x) = \frac{\psi_n(x)}{\zeta_n(x)} \ .

**Downward recursion for** :math:`A_n`. The upward recursion of :math:`A_n` is
unstable. LX-MIE therefore always uses the downward recursion

.. math::

   A_{n-1}(z) = \frac{n}{z} - \frac{1}{\dfrac{n}{z} + A_n(z)} \ ,

starting at :math:`n = N`. It is performed for both the complex argument
:math:`z = mx` and the real argument :math:`z = x`.

**Starting value from continued fractions.** The starting value
:math:`A_N(z)` is calculated with the method of continued fractions by Lentz
[Len76]_. It is based on the representation of the ratio of two Bessel
functions of consecutive order as a continued fraction,

.. math::

   A_N(z) = - \frac{N}{z} + \frac{J_{N-1/2}(z)}{J_{N+1/2}(z)}
          = - \frac{N}{z} + \left[ a_1; a_2, a_3, \ldots \right] \ ,
   \qquad
   a_k = \left( -1 \right)^{k+1} \frac{2 \left( N + k - 1/2 \right)}{z} \ ,

which Lentz's algorithm evaluates as a product of partial numerators and
denominators. The iteration is stopped when the relative difference between
the current partial numerator and denominator is smaller than
:math:`10^{-10}`. In contrast to the commonly used approach of starting the
downward recursion at an order much larger than :math:`N` with an arbitrary
value, this method does not require any additional terms and yields the
starting value with a known accuracy.

**Upward recursion for** :math:`B_n` **and** :math:`C_n`. These two
quantities are stable in upward direction and are calculated from

.. math::

   B_n(x) = - \frac{n}{x} + \frac{1}{\dfrac{n}{x} - B_{n-1}(x)} \ ,
   \qquad
   C_n(x) = C_{n-1}(x) \, \frac{B_n(x) + n/x}{A_n(x) + n/x} \ ,

with analytical starting values for :math:`B_0` and :math:`C_1` in terms of
trigonometric functions.

Angular functions
~~~~~~~~~~~~~~~~~

The angular functions :math:`\pi_n` and :math:`\tau_n` are calculated by the
upward recurrence relations of Wiscombe [Wis79]_. With :math:`\pi_0 = 0` and
:math:`\pi_1 = 1`,

.. math::

   s = \mu \, \pi_n \ , \qquad t = s - \pi_{n-1} \ ,

.. math::

   \pi_{n+1} = s + \frac{n + 1}{n} \, t \ , \qquad
   \tau_n = n \, t - \pi_{n-1} \ .

Legendre series
~~~~~~~~~~~~~~~

The coefficients of the :ref:`Legendre series <legendre-series>` of the phase
function are not obtained by a numerical integration over the scattering
angle, but directly from the Mie coefficients, following the procedure of
Dave [Dav70]_. The accuracy of the coefficients therefore does not depend
on the resolution of an angular quadrature.

Sign convention of the refractive index
---------------------------------------

The sign of the imaginary part of the refractive index depends on the chosen
time dependence of the electromagnetic wave. LX-MIE internally uses

.. math::

   m = n - \mathrm{i} k \ , \qquad k \geq 0 \ .

Every refractive index passed to LX-MIE is converted to this convention by
replacing the imaginary part with :math:`-|k|`. Consequently, either sign can
be used for the imaginary part. As a result of this convention, the
scattering amplitudes are the complex conjugates of those obtained with
:math:`m = n + \mathrm{i} k`, the convention used by Bohren & Huffman
[BH83]_.

Size distributions
------------------

The treatment of particle size distributions is described in
:doc:`size_distributions`.

References
----------

.. [KH18] Kitzmann, D. & Heng, K. (2018), *Optical properties of potential
   condensates in exoplanetary atmospheres*, MNRAS 475, 94–107,
   `arXiv:1710.04946 <https://arxiv.org/abs/1710.04946>`_

.. [BH83] Bohren, C. F. & Huffman, D. R. (1983), *Absorption and Scattering
   of Light by Small Particles*, Wiley

.. [vdH57] van de Hulst, H. C. (1957), *Light Scattering by Small
   Particles*, Wiley

.. [Mie08] Mie, G. (1908), *Beiträge zur Optik trüber Medien, speziell
   kolloidaler Metallösungen*, Annalen der Physik 330, 377–445

.. [Len76] Lentz, W. J. (1976), *Generating Bessel functions in Mie
   scattering calculations using continued fractions*, Applied Optics 15,
   668–671

.. [Wis79] Wiscombe, W. J. (1979), *Mie scattering calculations: Advances in
   technique and fast, vector-speed computer codes*, NCAR Technical Note
   NCAR/TN-140+STR

.. [Dav70] Dave, J. V. (1970), *Coefficients of the Legendre and Fourier
   series for the scattering functions of spherical particles*, Applied
   Optics 9, 1888–1896

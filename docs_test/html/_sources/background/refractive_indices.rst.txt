Refractive index compilation
============================

The ``compilation/`` directory of the repository contains the optical
constants of 32 potential condensates in planetary, exoplanetary, and other
astrophysical environments. The compilation and its data sources are described in
detail in Kitzmann & Heng (2018) [KH18]_.

.. note::

   The compilation is part of the repository, but not of the Python package
   on PyPI. To use it with ``pip install lxmie``, download the files from the
   `repository <https://github.com/newstrangeworlds/lx-mie/tree/master/compilation>`_.

File format
-----------

Each file contains three header lines and three tab-separated data columns:

- The first header line contains the chemical formula of the condensate and
  the references of the data, as well as additional information, e.g. on the
  temperature or on anisotropic and amorphous materials.
- The second header line contains the number of data points in the file.
- The third header line describes the three data columns.

The first data column is the wavelength in |micron|, the second one the real
part :math:`n` of the refractive index, and the third one the corresponding
imaginary part :math:`k`.

.. code-block:: text

   #Refractive index of amorphous (glassy) MgSiO3[s]; Dorschner et al. (1995), A&A Vol. 300; Database of Optical Constants for Cosmic Dust, Laboratory Astrophysics Group of the AIU Jena
   #Data points: 109
   #Wavelength (micron) n  k
   2.000000e-01	1.691660e+00	0.000000e+00
   2.200000e-01	1.659980e+00	0.000000e+00
   ...

The files can be used directly as the :ref:`refractive index file
<refractive-index-file>` of the command-line tool. Reading them with Python
is shown in :doc:`../python/examples`.

Available species
-----------------

The following table lists all files of the compilation with the wavelength
range that is covered by the data and the references given in the file
headers. If you use these data, please cite the original references in
addition to Kitzmann & Heng (2018).

.. csv-table::
   :file: ../_data/compilation.csv
   :header: "File", "Material", "Wavelengths (|micron|)", "Data points", "References"
   :widths: 20 18 12 8 42

Using other refractive indices
------------------------------

LX-MIE is not restricted to the species of the compilation. The Python
package and the C++ library accept arbitrary refractive indices. For the
command-line tool, any file that follows the format described in
:ref:`refractive-index-file` can be used.

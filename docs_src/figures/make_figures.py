"""
Creates the figures of the LX-MIE documentation.

Every figure is written twice, once for the light and once for the dark
colour scheme of the documentation theme.

Requirements: lxmie (>= 1.2), numpy, matplotlib

Usage:
    python make_figures.py <path to the LX-MIE repository>

The path to the LX-MIE repository is required for the refractive index
compilation. The figures are saved in the same directory as this script.
"""

import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import lxmie


output_dir = os.path.dirname(os.path.abspath(__file__))
repository = sys.argv[1] if len(sys.argv) > 1 else "."


themes = {
    "light": dict(
        surface="#fcfcfb", text="#0b0b0b", muted="#52514e",
        grid="#e1e0d9", axis="#c3c2b7",
        series=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark": dict(
        surface="#1a1a19", text="#f0efec", muted="#c3c2b7",
        grid="#2c2c2a", axis="#383835",
        series=["#3987e5", "#d95926", "#199e70"]),
}


def apply_theme(theme):
    plt.rcParams.update({
        "figure.facecolor": theme["surface"],
        "axes.facecolor": theme["surface"],
        "savefig.facecolor": theme["surface"],
        "text.color": theme["text"],
        "axes.labelcolor": theme["text"],
        "axes.edgecolor": theme["axis"],
        "xtick.color": theme["muted"],
        "ytick.color": theme["muted"],
        "grid.color": theme["grid"],
        "grid.linewidth": 0.6,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": plt.cycler(color=theme["series"]),
        "lines.linewidth": 2.0,
        "font.size": 10,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "legend.frameon": False,
        "figure.dpi": 100,
        "savefig.dpi": 200,
    })


def save(fig, name, theme_name):
    fig.savefig(os.path.join(output_dir, f"{name}_{theme_name}.png"))
    plt.close(fig)


def read_refractive_index(path):
    data = np.loadtxt(path, comments="#")
    data = np.atleast_2d(data)
    return list(data[:, 0]), list(data[:, 1]), list(data[:, 2])


# ---------------------------------------------------------------------------
# Figure 1: efficiencies as a function of the size parameter
# ---------------------------------------------------------------------------

def figure_efficiencies(theme_name):
    x = np.logspace(-1, 2, 3000)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5), sharey=True, constrained_layout=True)

    for ax, (n, k) in zip(axes, [(1.5, 0.01), (1.5, 0.1)]):
        results = lxmie.mie([n] * len(x), [k] * len(x), list(x))

        ax.plot(x, [r.q_ext for r in results], label=r"$Q_\mathrm{ext}$")
        ax.plot(x, [r.q_sca for r in results], label=r"$Q_\mathrm{sca}$")
        ax.plot(x, [r.q_abs for r in results], label=r"$Q_\mathrm{abs}$")

        ax.set_xscale("log")
        ax.set_xlabel(r"size parameter $x$")
        ax.set_title(f"m = {n} + {k}i")

    axes[0].set_ylabel("efficiency")
    axes[0].set_ylim(0, 4.6)
    axes[1].legend(loc="upper left")

    save(fig, "efficiencies", theme_name)


# ---------------------------------------------------------------------------
# Figure 2: phase functions for different size parameters
# ---------------------------------------------------------------------------

def figure_phase_function(theme_name):
    angles = np.linspace(0.0, 180.0, 721)
    cos_angles = list(np.cos(np.radians(angles)))

    fig, ax = plt.subplots(figsize=(6.5, 3.8), constrained_layout=True)

    for x in [0.1, 2.0, 10.0]:
        result, s1, s2 = lxmie.mie_amplitudes(1.33, 0.0, x, cos_angles)
        phase_function = lxmie.phase_function(s1, s2, x, result.q_sca)

        ax.plot(angles, phase_function, label=f"x = {x:g}")

    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_xticks(range(0, 181, 30))
    ax.set_xlabel("scattering angle (degrees)")
    ax.set_ylabel("phase function")
    ax.set_title("Phase function of a sphere with m = 1.33")
    ax.legend()

    save(fig, "phase_function", theme_name)


# ---------------------------------------------------------------------------
# Figure 3: the analytical size distributions
# ---------------------------------------------------------------------------

def figure_size_distributions(theme_name):
    distributions = [
        ("log-normal, $r_m = 1$, $\\sigma_g = 1.5$", lxmie.lognormal_distribution(1.0, 1.5)),
        ("gamma, $r_\\mathrm{eff} = 1$, $v_\\mathrm{eff} = 0.1$", lxmie.gamma_distribution(1.0, 0.1)),
        ("power law, $p = 3.5$, 0.05 – 5", lxmie.power_law_distribution(0.05, 5.0, 3.5)),
    ]

    fig, ax = plt.subplots(figsize=(6.5, 3.8), constrained_layout=True)

    for label, distribution in distributions:
        radius = np.array(distribution.radius)
        weight = np.array(distribution.weight)

        # the weights are n(r) dr on a grid equidistant in ln(r),
        # i.e. they are proportional to r n(r) = dN/dln(r)
        dn_dlnr = weight[1:-1]
        ax.plot(radius[1:-1], dn_dlnr / dn_dlnr.max(), label=label)

    ax.set_xscale("log")
    ax.set_xlim(0.03, 30)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("particle radius (micron)")
    ax.set_ylabel(r"$\mathrm{d}N / \mathrm{d}\ln r$ (scaled to maximum)")
    ax.set_title("Analytical size distributions")
    ax.legend()

    save(fig, "size_distributions", theme_name)


# ---------------------------------------------------------------------------
# Figure 4: single particle size compared to a size distribution
# ---------------------------------------------------------------------------

def figure_distribution_spectrum(theme_name):
    wavelengths, n_real, n_imag = read_refractive_index(
        os.path.join(repository, "compilation", "MgSiO3_amorph_glass.dat"))

    radius = 1.0
    size_parameters = [2.0 * math.pi * radius / wl for wl in wavelengths]
    single = lxmie.mie(n_real, n_imag, size_parameters)

    distribution = lxmie.lognormal_distribution(radius, 1.5)
    integrated = lxmie.mie_distribution(n_real, n_imag, wavelengths, distribution)

    fig, ax = plt.subplots(figsize=(6.5, 3.8), constrained_layout=True)

    ax.plot(wavelengths, [r.q_ext for r in single], label="single radius, 1 micron")
    ax.plot(wavelengths, [r.q_ext for r in integrated],
            label=r"log-normal, $r_m$ = 1 micron, $\sigma_g$ = 1.5")

    ax.set_xscale("log")
    ax.set_xlim(0.2, 100)
    ax.set_ylim(0, None)
    ax.set_xlabel("wavelength (micron)")
    ax.set_ylabel(r"extinction efficiency $Q_\mathrm{ext}$")
    ax.set_title(r"Amorphous MgSiO$_3$")
    ax.legend()

    save(fig, "distribution_spectrum", theme_name)


for theme_name, theme in themes.items():
    apply_theme(theme)

    figure_efficiencies(theme_name)
    figure_phase_function(theme_name)
    figure_size_distributions(theme_name)
    figure_distribution_spectrum(theme_name)

print("Figures written to", output_dir)

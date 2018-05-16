# Conquest_STM

This is a Python package to simulate Scanning Tunnelling Microscopy (STM) data under the Tersoff-Hamann approximation from the Conquest [1] code for Density Functional Theory (DFT) calculations.

## Motivation

STM is a vital tool in topographical and electronic studies of conducting and semi-conducting surfaces. The quantum tunnelling process that allows an STM to operate is not always straightforward. This often necessitates simulated images to aid in interpreting experimental data.

An important simulation technique is DFT - a computational method to calculate the electron density around atoms. Whilst a very successful method in simulations of condensed matter, it has some drawbacks, one of which is its scaling with system size, which goes as $\mathcal{O}(N^3)$. This limits systems of interest to roughly 1000 atoms even on high-performance computing systems.

This bottleneck is due to finding the Kohn-Sham eigenstates by diagonalising the Hamiltonian, or by variational minimisation of the energy and constraining the states to be orthogonal. By representing the wavefunction by pseudo-atomic orbitals (PAO) - orbitals localised around each atomic nucleus, rather than, say, plane waves - and solving for the density matrix, linear scaling is possible for DFT. Such a method is used in Conquest [1], a linear-scaling DFT code developed jointly by University College London (UCL) in the United Kingdom and the National Institute for Materials Science (NIMS) in Japan. This allows for systems of up to a million atoms to be accurately studied [2].

The Tersoff-Hamann approximation for calculating an STM current [3] does not explicitly include the atoms of the STM tip. Rather, only the wavefunction of the sample atoms at the position of the tip is required. This is not an issue for wavefunctions represented by plane waves, which can be accurate throughout the simulation cell. However, for PAOs which are only accurate around the sample atoms, this poses a problem for simulating STM data. A solution to this [4] is to evaluate the Bardeen integral [5] on a surface close to the sample atoms and propagate the wavefunctions up to the tip.

## Brief Mathematical Overview

Please note that atomic units are used for simplicity as well as consistency with CONQUEST.

As shown by Bardeen, the tunnelling current (what we want to simulate) is given by
$$
    I = 2\pi \sum_{\mu, \nu} \left| M_{\mu \nu} \right|^2 \delta(\varepsilon_{\mu} - \varepsilon_{\nu} - V),
$$

This sum is performed over the tip, $\mu$, and sample, $\nu$, states with energies $\varepsilon_i$, that lie in the bias voltage energy window, $V$. That is, $\varepsilon_F < \varepsilon_i < \varepsilon_F + V$ (or the appropriate inequality for a negative $V$) for the sample Fermi Level, $\varepsilon_F$.

The Bardeen integral is
$$
    M_{\mu \nu} = -\frac{1}{2} \iint_\Sigma \left[\chi_\mu^*(\mathbf{r}) \nabla\psi_\nu(\mathbf{r}) - \psi_\nu(\mathbf{r})\nabla\chi_\mu^*(\mathbf{r})\right] \cdot \mathop{}\!\mathrm{d}^2\mathbf{r},
$$

where $\chi$ and $\psi$ are the tip and sample wavefunctions respectively, and $\Sigma$ is ANY surface in the vacuum region between them. The method of Paz and Soler involves evaluating this integral close to the sample atoms, as it is $\psi(\mathbf{r})$ that becomes less accurate at larger distances. This wavefunction is calculated by CONQUEST during a simulation and is thus an input to this program (or rather the basis sets and coefficients are). The tip state must be determined by this program.

Explicit simulation of the atoms that constitute the STM probe tip is costly. Qualitative agreement with experiment can usually be achieved through the Tersoff-Hamann approximation. Representing the tip as a Dirac delta peak at $\mathbf{R}$ and solving the Schrödinger equation yields a Green's function that represents a spherical $s$-wave tip state:
$$
    \begin{aligned}
        \left(-\frac{1}{2} \nabla^2 - (\phi - \varepsilon_F) \right) G(\mathbf{r} - \mathbf{R}) & = -\delta(\mathbf{r} - \mathbf{R}), \\
        \left(\nabla^2 -  \kappa^2 \right) G(\mathbf{r} - \mathbf{R}) & = -2\delta(\mathbf{r} - \mathbf{R}); \\
        G(\mathbf{r} - \mathbf{R}) & = \frac{e^{-\kappa \left| \mathbf{r} - \mathbf{R} \right|}}{4 \pi \left|\mathbf{r} - \mathbf{R}\right|};
    \end{aligned}
$$

where the decay constant,
$$
    \kappa = \sqrt{2(\phi - \varepsilon_F)},
$$

is given by the work function of the sample, $\phi$.


Substitution of $\chi(\mathbf{r})$ with $G(\mathbf{r} - \mathbf{R})$ yields a simple expression for the Bardeen integral as a function of the tip position:
$$
    \begin{aligned}
        M_\varepsilon(\mathbf{R}) & = -\frac{1}{2} \iint_\Sigma \left[G^*(\mathbf{r} - \mathbf{R}) \nabla\psi_\varepsilon (\mathbf{r}) - \psi_\varepsilon(\mathbf{r})\nabla G^*(\mathbf{r} - \mathbf{R})\right] \cdot \mathop{}\!\mathrm{d}^2\mathbf{r}, \\
        & = - \sqrt{2 \pi \kappa} \psi_\varepsilon(\mathbf{R});
    \end{aligned}
$$

where $\psi_\varepsilon$ now represents the sample wavefunction of the state with energy $\varepsilon$.

Before $M_\varepsilon(\mathbf{R})$ can be substituted into our equation for $I$ we must take care of the Dirac delta in $I$ to get a computable expression. This delta function describes elastic tunnelling of the electrons from discrete states of exact equal energy. In a real system, thermal broadening of the states turns these discrete peaks into a continuous function. Furthermore, for a positive bias voltage, $V$, electrons will tunnel from the sample to the tip. We can thus replace the delta function with the Fermi-Dirac occupation factor, $f_{FD}(\varepsilon - V)$, where
$$
    f_{FD}(E) = \frac{1}{e^{(E - \varepsilon_F) kT} - 1},
$$

$k$ is Boltzmann's constant, and $T$ is the temperature of the sample. For the case of $V < 0$, it is unoccupied states in the sample that are tunnelled into by electrons in the tip, leading to the factor $[1 - f(\varepsilon - V)]$ instead. 

The tunnelling current at a given tip position is thus
$$
    \begin{aligned}
        I(\mathbf{R}) & = 2\pi \sum_{\varepsilon > \varepsilon_F}^{\varepsilon < \varepsilon_F + V} \left| M_\varepsilon(\mathbf{R}) \right|^2 f_{FD}(\varepsilon - V) & \quad \text{for} \quad V > 0, \\
        & = 2\pi \sum_{\varepsilon < \varepsilon_F}^{\varepsilon > \varepsilon_F + V} \left| M_\varepsilon(\mathbf{R}) \right|^2 [1 - f_{FD}(\varepsilon - V)] & \quad \text{for} \quad V < 0.
    \end{aligned}
$$

Thus, the tunnelling current may be obtained from the sample wavefunction, $\psi(\mathbf{r})$.

## Getting Started

An example script, ```example.py``` is included that shows the steps required from getting from Conquest output files to an STM image.

### Prerequisites

This project is written in Python 2.7. ```numpy```, ```matplotlib```, and ```skimage``` are required packages.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

1. <http://www.order-n.org>
2. D. R. Bowler and T. Miyazaki, J. Phys. Cond. Matter, **22**, 7, 074207 (2010)
3. J. Tersoff and D. R. Hamann, Phys. Rev. Lett., **50**, 1998-2001 (1983)
4. O. Paz and J. M. Soler, Phys. Stat. Sol. B, **5**, 1080-1094 (2005)
5. J. Bardeen, Phys. Rev. Lett., **6**, 57-59 (1961)

## Acknowledgments

* Prof. David Bowler, University College London, London Centre for Nanotechnology, for supervising this Master's project.

The $\text{\LaTeX}$ equations in this document were compile with ```readme2tex``` from the repository <https://github.com/leegao/readme2tex>.

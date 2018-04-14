# Conquest_STM

This is a Python package to simulate Scanning Tunnelling Microscopy (STM) data under the Tersoff-Hamann approximation from the Conquest [1] code for Density Functional Theory (DFT) calculations.

## Motivation

STM is a vital tool in topographical and electronic studies of conducting and semi-conducting surfaces. The quantum tunnelling process that allows an STM to operate is not always straightforward. This often necessitates simulated images to aid in interpreting experimental data.

An important simulation technique is DFT - a computational method to calculate the electron density around atoms. Whilst a very successful method in simulations of condensed matter, it has some drawbacks, one of which is its scaling with system size, which goes as O(N<sup>3</sup>). This limits systems of interest to roughly 1000 atoms even on high-performance computing systems.

This bottleneck is due to finding the Kohn-Sham eigenstates by diagonalising the Hamiltonian, or by variational minimisation of the energy and constraining the states to be orthogonal. By representing the wavefunction by pseudo-atomic orbitals (PAO) - orbitals localised around each atomic nucleus, rather than, say, plane waves - and solving for the density matrix, linear scaling is possible for DFT. Such a method is used in Conquest [1], a linear-scaling DFT code developed jointly by University College London (UCL) in the United Kingdom and the National Institute for Materials Science (NIMS) in Japan.

The Tersoff-Hamann approximation for calculating an STM current does not explicitly include the atoms of the STM tip. Rather, only the wavefunction of the sample atoms at the position of the tip is required. This is not an issue for wavefunctions represented by plane waves, which can be accurate throughout the simulation cell. However, for PAOs which are only accurate around the sample atoms, this poses a problem for simulating STM data. A solution to this is to evaluate the Bardeen integral on a surface close to the sample atoms and propagate the wavefunctions up to the tip.

## Getting Started

An example script, ```example.py``` is included that shows the steps required from getting from Conquest output files to an STM image.

### Prerequisites

This project is written in Python 2.7. ```numpy```, ```matplotlib```, and ```skimage``` are required packages.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

1. <http://www.order-n.org>
2. D. R. Bowler and T. Miyazaki, J. Phys. Cond. Matter, **22**, 7, 074207 (2010)
3. O. Paz and J. M. Soler, Phys. Stat. Sol. B, **5**, 1080-1094 (2005)

## Acknowledgments

* Prof. David Bowler, University College London, London Centre for Nanotechnology, for supervising this Master's project.

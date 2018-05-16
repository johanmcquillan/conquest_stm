# Conquest_STM

This is a Python package to simulate Scanning Tunnelling Microscopy (STM) data under the Tersoff-Hamann approximation from the Conquest [1] code for Density Functional Theory (DFT) calculations.

## Motivation

STM is a vital tool in topographical and electronic studies of conducting and semi-conducting surfaces. The quantum tunnelling process that allows an STM to operate is not always straightforward. This often necessitates simulated images to aid in interpreting experimental data.

An important simulation technique is DFT - a computational method to calculate the electron density around atoms. Whilst a very successful method in simulations of condensed matter, it has some drawbacks, one of which is its scaling with system size, which goes as <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/3c5638d37c66029a01e04817cbad1d37.svg?invert_in_darkmode" align=middle width=48.557355pt height=26.70657pt/>. This limits systems of interest to roughly 1000 atoms even on high-performance computing systems.

This bottleneck is due to finding the Kohn-Sham eigenstates by diagonalising the Hamiltonian, or by variational minimisation of the energy and constraining the states to be orthogonal. By representing the wavefunction by pseudo-atomic orbitals (PAO) - orbitals localised around each atomic nucleus, rather than, say, plane waves - and solving for the density matrix, linear scaling is possible for DFT. Such a method is used in Conquest [1], a linear-scaling DFT code developed jointly by University College London (UCL) in the United Kingdom and the National Institute for Materials Science (NIMS) in Japan. This allows for systems of up to a million atoms to be accurately studied [2].

The Tersoff-Hamann approximation for calculating an STM current [3] does not explicitly include the atoms of the STM tip. Rather, only the wavefunction of the sample atoms at the position of the tip is required. This is not an issue for wavefunctions represented by plane waves, which can be accurate throughout the simulation cell. However, for PAOs which are only accurate around the sample atoms, this poses a problem for simulating STM data. A solution to this [4] is to evaluate the Bardeen integral [5] on a surface close to the sample atoms and propagate the wavefunctions up to the tip.

## Brief Mathematical Overview

Please note that atomic units are used for simplicity as well as consistency with CONQUEST.

As shown by Bardeen, the tunnelling current (what we want to simulate) is given by
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/2c5bc2bdbc18e99f66c7157c157089f9.svg?invert_in_darkmode" align=middle width=239.811pt height=38.38758pt/></p>

This sum is performed over the tip, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode" align=middle width=9.86799pt height=14.10255pt/>, and sample, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/b49211c7e49541e500c32b4d56d354dc.svg?invert_in_darkmode" align=middle width=9.132585pt height=14.10255pt/>, states with energies <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/e76eead5067c5c0fdf70614a08ab0c95.svg?invert_in_darkmode" align=middle width=12.27039pt height=14.10255pt/>, that lie in the bias voltage energy window, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/a9a3a4a202d80326bda413b5562d5cd1.svg?invert_in_darkmode" align=middle width=13.192575pt height=22.38192pt/>. That is, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/8a5f6e4f652372d9bbf551ffd6d7f8dd.svg?invert_in_darkmode" align=middle width=127.215495pt height=22.38192pt/> (or the appropriate inequality for a negative <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/a9a3a4a202d80326bda413b5562d5cd1.svg?invert_in_darkmode" align=middle width=13.192575pt height=22.38192pt/>) for the sample Fermi Level, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/391eb7f9ef426346bab88aabcda164f6.svg?invert_in_darkmode" align=middle width=17.705325pt height=14.10255pt/>.

The Bardeen integral is
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/941e8651fc1f6369950f358399533a79.svg?invert_in_darkmode" align=middle width=374.68695pt height=37.35204pt/></p>

where <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/c91091e68f0e0113ff161179172813ac.svg?invert_in_darkmode" align=middle width=10.246995pt height=14.10255pt/> and <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/7e3c241c2dec821bd6c6fbd314fe4762.svg?invert_in_darkmode" align=middle width=11.255475pt height=22.74591pt/> are the tip and sample wavefunctions respectively, and <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode" align=middle width=11.82786pt height=22.38192pt/> is ANY surface in the vacuum region between them. The method of Paz and Soler involves evaluating this integral close to the sample atoms, as it is <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/2de90df5ef3c9406edcd8e1d78a9e5e7.svg?invert_in_darkmode" align=middle width=31.75161pt height=24.56553pt/> that becomes less accurate at larger distances. This wavefunction is calculated by CONQUEST during a simulation and is thus an input to this program (or rather the basis sets and coefficients are). The tip state must be determined by this program.

Explicit simulation of the atoms that constitute the STM probe tip is costly. Qualitative agreement with experiment can usually be achieved through the Tersoff-Hamann approximation. Representing the tip as a Dirac delta peak at <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.125155pt height=22.473pt/> and solving the Schrödinger equation yields a Green's function that represents a spherical <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode" align=middle width=7.6767405pt height=14.10255pt/>-wave tip state:
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/9bfdbaa8c533b512ad4107ffc190d94b.svg?invert_in_darkmode" align=middle width=331.9932pt height=113.807595pt/></p>

where the decay constant,
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/077a05818b3243b70fe48a1b3b31feec.svg?invert_in_darkmode" align=middle width=121.777095pt height=19.654965pt/></p>

is given by the work function of the sample, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode" align=middle width=9.757935pt height=22.74591pt/>.


Substitution of <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/84fa817903364a747cabef21126476d0.svg?invert_in_darkmode" align=middle width=30.740985pt height=24.56553pt/> with <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/7bd3ed61c120134b1f6edf60f44cd32c.svg?invert_in_darkmode" align=middle width=67.53879pt height=24.56553pt/> yields a simple expression for the Bardeen integral as a function of the tip position:
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/565061f3a2762e9d87da52e3531160b3.svg?invert_in_darkmode" align=middle width=457.01205pt height=63.73719pt/></p>

where <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/72f2218c93bfb88863a87d988cfcdb1a.svg?invert_in_darkmode" align=middle width=16.855245pt height=22.74591pt/> now represents the sample wavefunction of the state with energy <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/9ae7733dac2b7b4470696ed36239b676.svg?invert_in_darkmode" align=middle width=7.6369095pt height=14.10255pt/>.

Before <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/ee0de26e7f8622fc36aee8ce81bfbfaf.svg?invert_in_darkmode" align=middle width=49.783305pt height=24.56553pt/> can be substituted into our equation for <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/21fd4e8eecd6bdf1a4d3d6bd1fb8d733.svg?invert_in_darkmode" align=middle width=8.4843pt height=22.38192pt/> we must take care of the Dirac delta in <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/21fd4e8eecd6bdf1a4d3d6bd1fb8d733.svg?invert_in_darkmode" align=middle width=8.4843pt height=22.38192pt/> to get a computable expression. This delta function describes elastic tunnelling of the electrons from discrete states of exact equal energy. In a real system, thermal broadening of the states turns these discrete peaks into a continuous function. Furthermore, for a positive bias voltage, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/a9a3a4a202d80326bda413b5562d5cd1.svg?invert_in_darkmode" align=middle width=13.192575pt height=22.38192pt/>, electrons will tunnel from the sample to the tip. We can thus replace the delta function with the Fermi-Dirac occupation factor, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/c63e32227981794f696b8263edb42426.svg?invert_in_darkmode" align=middle width=83.672325pt height=24.56553pt/>, where
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/221a083892ea5917e0c9f4bb8c4f2bee.svg?invert_in_darkmode" align=middle width=186.1926pt height=34.61007pt/></p>

<img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.041505pt height=22.74591pt/> is Boltzmann's constant, and <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode" align=middle width=11.84502pt height=22.38192pt/> is the temperature of the sample. For the case of <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/2b896cf26468d63c705295ac622e1c08.svg?invert_in_darkmode" align=middle width=43.26465pt height=22.38192pt/>, it is unoccupied states in the sample that are tunnelled into by electrons in the tip, leading to the factor <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/240808c6aadbbbfad2108d0320585b2c.svg?invert_in_darkmode" align=middle width=100.74174pt height=24.56553pt/> instead. 

The tunnelling current at a given tip position is thus
<p align="center"><img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/d777fedd83e7062a62423228d3331c39.svg?invert_in_darkmode" align=middle width=441.1176pt height=108.15387pt/></p>

Thus, the tunnelling current may be obtained from the sample wavefunction, <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/2de90df5ef3c9406edcd8e1d78a9e5e7.svg?invert_in_darkmode" align=middle width=31.75161pt height=24.56553pt/>.

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

The <img src="https://rawgit.com/johanmcquillan/conquest_stm/master/svgs/c068b57af6b6fa949824f73dcb828783.svg?invert_in_darkmode" align=middle width=42.05817pt height=22.407pt/> equations in this document were compile with ```readme2tex``` from the repository <https://github.com/leegao/readme2tex>.

---
title: 'CaNS-Fizzy: A GPU-accelerated finite difference solver for turbulent two-phase flows'
tags:
  - Direct Numerical Simulation
  - Multiphase Flow
  - GPU
authors:
  - name: Giandomenico Lupo
    corresponding: true
    orcid: 0000-0002-1095-118X
    equal-contrib: true
    affiliation: 1
  - name: Peter Wellens
    orcid: 0000-0003-0009-9058
    affiliation: 2
  - name: Pedro Costa
    orcid: 0000-0001-7010-1040
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Delft University of Technology, Department of Process \& Energy, Delft, The Netherlands
   index: 1
 - name: Delft University of Technology, Department of Marine \& Transport Technology, Delft, The Netherlands
   index: 2
date: 30 October 2024
bibliography: paper.bib

---

# Summary

CaNS-Fizzy -- Fizzy for short -- is a GPU-accelerated numerical solver for massively-parallel Direct Numerical Simulations (DNS) of incompressible two-phase flows. A DNS enables direct access to all flow quantities, resolved in time and space at all relevant continuum scales. The resulting numerical experiments provide complete data sets for the analysis of the detailed mechanisms underlying the flow, particularly the interaction between the chaotic and multi-scale dynamics of turbulence and the interface movement and deformation. The insights gained can guide the design and operation of various applications, such as boiling heat transfer, liquid-liquid extraction, gas-liquid reactors, absorption and stripping columns, distillation columns, liquid combustion appliances, in all of which the rate of heat and mass transfer between phases is proportional to the interfacial area. Fizzy's two-phase capabilities were implemented using the efficient, GPU-accelerated Navier-Stokes solver CaNS [@Costa:2018] as a base.

# Statement of need

Fizzy is suited for large-scale direct numerical simulations of canonical incompressible two-phase flows, from simple laminar cases to bubble/droplet-laden turbulent suspensions. These flows may be computationally expensive due to the stringent resolution requirements imposed by the direct solution of immersed interfaces dispersed throughout the domain. This demands efficient use of the capabilities of modern computing systems; Fizzy has been developed to include key desirable features that enable this objective: a one-fluid formulation of the two-phase flow governing equations, the use of a fast direct solver for the pressure Poisson equation, and an efficient distributed GPU porting with an interface capturing strategy that is suitable for GPU acceleration.
In addition to the momentum transfer and interface capturing, the code has the capability to solve heat transfer in both fluid phases, and thermal convection based on the Oberbeck-Boussinesq approximation.
Finally, the code has been extensively validated with several benchmark cases that demonstrate the different features of the solver, which are incorporated in the continuous integration workflows of the repository.
The GPU capabilities differentiate Fizzy from other commonly used open-source state-of-the-art incompressible two-phase flow solvers such as [boilingfoam](https://github.com/fmuni/boilingFoam-PUBLIC) [@boilingfoam] and [Basilisk](http://basilisk.fr/) [@gerris,@basilisk], which are more suited to smaller-scale direct numerical simulations in complex geometries. The recently published [FLOW36](https://github.com/MultiphaseFlowLab/FLOW36) [@flow36] is another GPU-ready code with similar features to Fizzy, differing in the interface capturing scheme and in the use of a pseudo-spectral instead of a finite difference approach.

# Mathematical model

A one-fluid formulation of the two-phase flow is employed, solving a single set of governing equations for both phases in the whole domain. The incompressible Navier-Stokes equation, the heat transport equation and the phase indicator transport equation are evolved in time to compute the velocity and pressure, temperature and phase indicator fields respectively. The latter identifies the regions of the domain occupied by either phase: it is continuous and smooth over the whole domain, and the interface between phases is diffuse. The thermophysical and transport properties (density, viscosity, thermal conductivity, specific heat capacity) are linearly mapped over the phase indicator field, and thus also continuous and smooth. The surface tension at the interface is included as a Continuous Surface Force (CSF) [@Brackbill:1992] in the Navier-Stokes equation.

# Methods and Implementation strategy

The governing equations are spatially discretized with a second-order finite difference scheme on a 3D Cartesian grid; a staggered grid arrangement is used for the velocity field, while all other quantities are stored at the cell centers; time integration is based on a low-storage three-step Runge-Kutta scheme. The incompressible Navier-Stokes equation is solved with a pressure correction scheme to enforce mass conservation, which yields a variable coefficient Poisson equation for the pressure correction: a splitting technique adapted from [@Dong:2012] and [@Dodd:2014] transforms this equation into a constant coefficient Poisson equation [@Frantzis:2019], enabling the use of the fast direct FFT solver of the CaNS code. Fizzy also allows for solving the conventional variable-coefficients problem using a geometric multigrid method through the [Hypre](https://github.com/hypre-space/hypre) library.
The interface capturing scheme for the phase indicator transport equation can be chosen between the Accurate Conservative Diffuse Interface (ACDI) scheme [@Jain:2022] and a tailored flavour of the THINC algebraic Volume-of-Fluid method [@thinc]. Both methods share a diffuse interface representation of the phase interface that allows for continuous and smooth mapping of the physical fields across the interface. The ACDI method requires no explicit interface reconstruction thanks to an interface regularization flux; for the VoF method the interface geometry in each cell is simplified to allow analytic calculation of the interface-cell intersection and of the advection fluxes at cell faces: thus in both methods the computational load is kept constant regardless of the local interface topology, making the algorithm particularly suited for parallelization on GPU architecture. For the ACDI method, the momentum equation includes the flux associated with the diffuse interface regularization, which allows for a mass--momentum consistent discretization and enables stability at high density contrasts between phases. The heat equation similarly includes the enthalpy flux associated with the interface regularization.

The code is written in modern Fortran, and is parallelized using MPI and OpenACC directives for GPU kernel offloading and host/device data transfer. As in CaNS, Fizzy leverages [cuDecomp](https://github.com/NVIDIA/cuDecomp) [@Romero:2022] for distributed memory calculations in pencil domain decompositions, and [cuFFT](https://docs.nvidia.com/cuda/cufft/) for computing Fourier transforms. These libraries are designed to work on NVIDIA GPUs; in the future Fizzy will support other GPU hardware, following updates on CaNS. On CPUs, the code uses [2DECOMP&FFT](https://github.com/2decomp-fft/2decomp-fft) [@Rolfo:2023] and [FFTW](https://www.fftw.org/) to perform the same operations.

Users can design and run a simulation by specifying the physical and computational parameters in a simple Fortran namelist input file. The code uses a modular, procedural design which makes extensions with different numerical methods or physical phenomena easy to develop. In the short term, we aim to allow for alternative schemes for spatial and temporal discretization, and introduce different interface capturing schemes, e.g. geometric VoF.

Finally, the code was designed so that important new computational features in the parent solver CaNS (e.g. porting efforts to other architectures) are easily propagated to Fizzy.

# Examples

\autoref{fig:examples} illustrates examples of two-phase flows simulated with this solver. The left panel shows how forced homogeneous isotropic turbulence breaks and deforms the interface of a liquid-liquid emulsion in a triperiodic domain; the right panel shows a hot gas bubble rising in a cold liquid, attaining the typical skirt shape while cooling down and simultaneously heating up the surrounding liquid in its wake.

![(Left) Simulation of a liquid-liquid emulsion in a three-periodic domain with sustained homogeneous isotropic turbulence. The colour represents the velocity magnitude. A diagonal plane with in-plane velocity vectors is shown. (Right) Three successive snapshots of a hot gas bubble rising in a cold liquid. A coloured volume rendering of the temperature field is shown both in the gas bubble and in the liquid wake. \label{fig:examples}](examples.png){ width=100% }

# Computational performance

Fizzy is tailored for large-scale simulations that exploit the computational capacity of modern GPU clusters with full GPU occupancy. The most relevant metric of the parallel efficiency for such scenario is a weak scaling test that determines the penalty in increased wall-clock time occurring when the problem size is increased alongside the computational resources. The liquid-liquid emulsion in homogeneous isotropic turbulence case of \autoref{fig:examples} (Left) has been used for this test: the size of the computational domain has been extended in one direction linearly with the number of GPU nodes employed. The test has been carried out on the GPU partition of the supercomputer Leonardo from Cineca, Italy; each computing node is equipped with four NVIDIA A100 SXM6 64GB GPUs, and is able to fit at full memory a $1024^3$ ($\sim$ 1 billion grid cells) computational box. \autoref{fig:performance} shows the performance penalty for the ACDI method, as the problem domain size (i.e. the number of spatial degrees of freedom) is increased from occupying 4 nodes (16 GPUs) to 64 nodes (256 GPUs): the 16 times larger computation takes only about 1.7 times longer than the original 4-node computation. A similar performance is obtained using the algebraic VoF method. The key contributor to the parallel performance is the interface capturing approach which prevents thread divergence in GPU kernels, as the computational load is independent of the local interface morphology. Indeed, very little sensitivity of the wall-clock time per iteration to the amount of interface area is observed, even for unsteady evolution of the interface with drastic topology changes during break-up events.

![Weak scaling performance on GPU nodes at full memory. The vertical axis shows the wall-clock time normalized by the four-node case.\label{fig:performance}](weak_scaling.png){ width=60% }

# Acknowledgements

We would like to thank Jordi Poblador-Ibanez, Suhas Jain, Naoki Hori, and Bergmann Óli Aðalsteinsson for insightful discussions that led to practical improvements of the numerical method.
This work was partially supported by a Cohesion Grant from the Mechanical Engineering Faculty at TU Delft awarded to Pedro Costa and Peter Wellens.

# References

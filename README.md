# R-Soft-Inverted-Pendulum
Modeling and control of the Soft Inverted Pendulum system with an actuated revolute joint at the base.

## Contents:
* [1. Requirements](#1-requirements)
* [2. Code Guide](#2-code-guide)
* [3. Modeling](#3-modeling)
* [4. Structural Properties](#4-structural-properties)
* [5. Control](#5-control)
* [6. Adaptive Control](#6-adaptive-control)
* [7. Simulations](#7-simulations)

## 1) Requirements
To run the scripts, you need to download Peter Corke's Robotics Toolbox, available at the following [link](https://petercorke.com/toolboxes/robotics-toolbox/). Additionally, to run the Simulink files, you need to have MATLAB version R2021b.

## 2) Code Guide
The package is organized into three directories that contain functions useful for the analysis and control of the model. Below is a description of these directories:

- `Della Santina`: This directory contains the functions provided with the project. These functions are used to compare the original model from the paper with the MATLAB implementation.

- `my_functions`: This directory contains math utility functions.

- `origin_soft_pendulum`: This directory contains MATLAB-generated functions that implement the dynamics of the Soft Inverted Pendulum.

The remaining scripts and functions are described in the following sections.

## 3) Modeling
- The modeling of the Soft Inverted Pendulum system is performed by the script `Soft_origin.m`. After calculating the relevant dynamic elements, it generates functions useful for simulations.

- The modeling of the R-Soft Inverted Pendulum system is implemented by the script `Rsoft_model.m`. Similarly to the previous script, it calculates the relevant dynamic elements and generates functions for simulations. Additionally, the script performs a linearization of the system.

- The equilibrium analysis is performed by the script `equilibria_analysis.m`, which analyzes the stability, reachability, and observability of the linearized system. The equilibria are derived from the `.mat` files generated by the `Rsoft_model.m` script.

- The scripts `validazione_coriolis.m` and `validazione_dinamica.m` compare the values of the dynamic matrices generated by the previous scripts, applying them to a set of state variable values.

- A purely kinematic simulation is implemented by the script `test_kin.m`.

## 4) Structural Properties
For a more in-depth control of the system's structural properties, an iterative calculation of the accessibility distribution and observability codistribution is performed. These distributions are calculated in the script `acc_obs_dist.m`. Additionally, several useful functions have been implemented:
- `lieBracket.m`: Implements the calculation of the Lie Bracket between vectors.
- `rowLieBracket.m`: Implements the calculation of the Lie Bracket between a vector and a covector.
- `filtration.m`: Implements the iterative calculation of filtration, useful for the accessibility distribution.
- `rowFiltration.m`: Implements the iterative calculation of filtration for the observability codistribution.

## 5) Control
For control purposes, the script `feedback_lin.m` has been implemented to calculate the inputs and new state, setting the desired output to be controlled. Additionally, the functions `collocatedFL.m` and `collocatedFL2.m` implement feedback linearization on \( \theta_r \) and \( \alpha_s \), respectively. The controller implementation can be found in the Simulink file `R_soft_sim.slx`, which has a dedicated section for it.

## 6) Adaptive Control
Adaptive control is implemented by the function `regressorSoftInverted.m`. This function returns the regressor of the Soft Inverted Pendulum as indicated in this [article](https://ieeexplore.ieee.org/abstract/document/9482817).

## 7) Simulations
The Simulink file `R_soft_sim.slx` contains several simulations:
- **Collocated Feedback Linearization Spong Like**: Control on the base joint variable.

- **Feedback Linearization Alpha**: Control on the tip inclination with respect to the fixed reference frame.

- **Autonomous Soft Inverted Pendulum**: Autonomous Soft Inverted Pendulum system or system with slowly varying inputs.

- **Adaptive Controller Soft Inverted Pendulum**: Adaptive control of the Soft Inverted Pendulum system. The results are fully comparable with those of the [article](https://ieeexplore.ieee.org/abstract/document/9482817).

The script `test_din.m` declares the simulation parameters and records the simulation values. The `plot_Rsoft.m` function has been implemented to visualize the system and its trajectory.


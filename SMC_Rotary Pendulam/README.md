# SMC Rotary Pendulum Control

## Overview

This repository contains the implementation of Sliding Mode Control (SMC) for a rotary inverted pendulum system. The project includes modeling, controller design, and simulation of the rotary pendulum dynamics using MATLAB and Simulink.

## Contents

- `smc_rotary_pendulam_setup.m`: MATLAB script that sets up parameters and runs simulations of the rotary pendulum under the designed Sliding Mode Controller.
- `RotaryPendulum.slx`: Simulink model of the rotary inverted pendulum system with sliding mode control (not uploaded).
- `Report.pdf`: Detailed report describing the system modeling, control design methodology, implementation details, and simulation results.

## Description

The code implements a robust sliding mode controller to stabilize the inherently unstable rotary inverted pendulum. It uses nonlinear system dynamics, accounts for model uncertainties, and applies SMC to ensure convergence despite disturbances.

The Simulink model captures the system dynamics and controller in a block-diagram form for simulation and visualization.

## Usage

- Run the MATLAB script `smc_rotary_pendulam_setup.m` to simulate the Rotary Pendulum system with SMC and generate performance plots.
- Open and simulate the `RotaryPendulum.slx` Simulink model to visualize and experiment with the system in a graphical environment.

---


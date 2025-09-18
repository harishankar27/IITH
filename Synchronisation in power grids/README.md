# Synchronisation in Power Grids


This project analyzes power grid synchronization using a generalized Kuramoto model, which models generators and loads as coupled nonlinear oscillators. Maintaining synchronization in power grids is critical to prevent failures and blackouts. The Kuramoto model provides a theoretical framework to study phase dynamics and grid stability.

- Simulates a simplified 3-node grid with two generators and one machine.
- Studies effects of perturbations in power demand on synchronization and power delivery.

-`network_dynamics.ipynb`: The notebook includes code for simulating oscillator phases, computing order parameters, and visualizing transient and steady-state network behavior.

## Results and Observations

- Under normal operation and moderate perturbations, the system returns to synchronized phase locking.
- Near critical perturbation, oscillatory and unstable synchronization occurs.
- Above threshold perturbations cause loss of synchronization with erratic power delivery.


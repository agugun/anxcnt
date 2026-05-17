# AXSCNT Documentation

Welcome to the AXSCNT technical documentation. This directory contains detailed architectural designs and module-specific diagrams.

## Core Architecture

- [Global System Architecture](architecture.md): The high-level contracts and engine orchestration.

## Physics Modules

Each module implements the core [Simulation Contracts](architecture.md#top) for specific physical domains.

### Classical Mechanics
- [Harmonic Oscillator](modules/oscillator.md)

### Thermodynamics
- [1D Heat Conduction (Implicit)](modules/heat_1d.md)
- [2D Heat Conduction (Implicit)](modules/heat_2d.md)

### Fluid Dynamics
- [Burgers Equation](modules/burgers.md)
- [Incompressible Flow (FEM)](modules/fluid_dynamics.md)

### Wave Phenomena
- [1D Wave Equation](modules/wave_1d.md)
- [2D Wave Equation](modules/wave_2d.md)

### Petroleum Engineering
- [1D Pressure Diffusivity](modules/pressure_1d.md)
- [Black Oil (2D/3D)](modules/black_oil.md)
- [Material Balance (MBA)](modules/mba.md)

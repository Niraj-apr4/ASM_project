# Cantilever Beam Simulation with Magnetic Actuation

A numerical simulation of **field-responsive metamaterials** — structures that deform under magnetic fields. Built using Discrete Differential Geometry (DDG) principles and implicit Euler time integration.

> 🚧 **Status:** 1D cantilever beam simulation complete. 2D grid extension in progress.

## Overview

This project models a cantilever beam as a chain of point masses connected by springs, subjected to gravity and an external magnetic force. The goal is to study how magnetic field direction and intensity, spring stiffness, and damping affect beam deflection — a step toward engineering adaptive, magnetically tunable structures.

## 1D Implementation

- **Discretization:** 11 lumped point masses connected by 10 linear spring-damper elements
- **Beam length:** 1.0 m, clamped at one end, free at the other
- **Forces modeled:**
  - Spring restoration (Hooke's law)
  - Gravity
  - Magnetic actuation (applied at the free end)
- **Time integration:** Implicit Euler method
- **Damping:** Simple viscous damping to dissipate energy each timestep

### Parametric Studies

| Study | What was varied |
|---|---|
| Directional loading | Magnetic field direction (NE, NW, SE, SW) |
| Field intensity | Magnetic susceptibility χ (10–100) |
| Spring stiffness | k from 10³ to 10⁵ N/m |
| Damping | Damping coefficient from 0.5 to 1.0 |

**Key result:** Deflection scales linearly with magnetic field intensity and inversely with spring stiffness, confirming the model behaves physically as expected.

## 2D Implementation — In Progress 🚧

Extending the framework to a full 11×11 node grid (121 nodes, 200 spring connections) to capture planar deformation patterns under combined gravitational and magnetic loading. Early results are promising — full parametric analysis and validation still underway.

## Simulation Parameters (1D)

| Parameter | Value |
|---|---|
| Spring stiffness (k) | 10⁴ N/m |
| Damping coefficient | 0.9 |
| Magnetic susceptibility (χ) | 50 |
| Node mass | 1.0 kg |
| Timestep | 0.001 s |
| Total simulation time | 1.0 s |

## Authors

- Niraj Kr Singh (AM25M807)
- V. Lalu (AM25M008)

## Future Work

- Complete 2D validation and parametric sweeps
- Extend to 3D magnetoelastic structures
- Nonlinear material models
- Experimental validation

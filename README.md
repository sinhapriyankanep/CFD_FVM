# 🔥 1D Heat Conduction in a Rod (Transient Analysis)

This project solves the **one-dimensional transient heat conduction problem** in a rod using the heat equation with Dirichlet boundary conditions and a non-linear initial temperature profile.



## 📘 Problem Statement

Consider a 1D rod of normalized length \( L = 1 \) m, governed by the transient heat equation:

\[
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}
\]

Where:
- \( T(x, t) \): Temperature distribution in space and time
- \( \alpha = 1 \): Thermal diffusivity (normalized)



### 🔒 Boundary Conditions

\[
T(0, t) = T_0, \quad T(1, t) = T_L
\]

---

## ⏳ Initial Condition

\[
T(x, 0) = T_0 + (T_L - T_0) \left( \frac{x}{L} \right)^2, \quad x \in (0, L)
\]

---

## 💻 Numerical Solution (Finite Difference Method)

This repository includes a numerical solution using the finite difference method with:

- Explicit or implicit time integration
- Fixed time step and spatial grid size
- Stability check (for explicit methods: CFL condition)

---

## 📁 Files

| File | Description |
|------|-------------|
|`1D_unsteady_diffusion_v2.m` | Code to solve the PDE using FVM |
| `README.md` | Problem description and repository guide |
| `CFD_A1` | PDF with the formulation and solution |




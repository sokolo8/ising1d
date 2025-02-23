---
layout: default
title: Ising Model Simulation
---

# 1D Ising Model with Transverse Field

The Hamiltonian of the system is:

$$ H = - \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - g \sum_{i=1}^{N} \sigma^x_i $$

where:
- \( N \) is the number of spin-1/2 particles.
- \( g \) is the transverse magnetic field strength.
- \( \sigma^x_i, \sigma^z_i \) are Pauli matrices.

## Installation
Run:
```bash
pip install numpy scipy matplotlib

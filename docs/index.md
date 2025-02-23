---
layout: default
title: 1D Ising model simulations
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.js">
</script>

# 1D Ising Model with Transverse Field

The Hamiltonian of the system is:
\begin{align}
H = - \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - g \sum_{i=1}^{N} \sigma^x_i
\end{align}
where:
- $N$ is the number of spin-1/2 particles.
- $g$ is the transverse magnetic field strength.
- $\sigma^x_i$, $\sigma^z_i$ spin-1/2 operators associated with site $i$ in the direction of $\hat{x}$, $\hat{z}$ respectively

## Installation
Run:
```bash
pip install numpy scipy matplotlib

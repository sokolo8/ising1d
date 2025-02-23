---
layout: default
title: 1D Ising model simulations
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.js">
</script>

# 1D Ising Model with Transverse Field

The Hamiltonian of the system is:

$$H = - \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - g \sum_{i=1}^{N} \sigma^x_i$$

where:
- $$N$$ is the number of spin-1/2 particles.
- $$g$$ is the transverse magnetic field strength.
- $$\sigma^x_i$$, $$\sigma^z_i$$ spin-1/2 operators associated with site $$i$$ in the direction of $$\hat{x}$$, $$\hat{z}$$ respectively

The model is in a $\mathbf{paramagnetic\ phase}$, when $$g>>1$$. In this case ground state of the system has no degeneracy and it is a state in which all spins are pointing along $$\hat{x}$$:

$$\ket{GS}=\ket{\rightarrow\rightarrow\rightarrow...\rightarrow\rightarrow\rightarrow}$$

On the other hand, the most excited state in this phase would be the following one:

$$\ket{\leftarrow\leftarrow\leftarrow...\leftarrow\leftarrow\leftarrow}$$

The model is in a \textbf{ferromagnetic phase}, when $$g<<1$$ or $$g \to 0$$. A ground state of the system has a two-fold degeneracy and it can be one of the following states (or their superposition):

$$\ket{\uparrow\uparrow\uparrow...\uparrow\uparrow\uparrow}, \quad \quad \ket{\downarrow\downarrow\downarrow...\downarrow\downarrow\downarrow}$$

On the other hand, the most excited state in this phase would be:

$$\ket{\uparrow\downarrow\uparrow...\downarrow\uparrow\downarrow} \quad \mathrm{or} \quad \ket{\downarrow\uparrow\downarrow...\uparrow\downarrow\uparrow}$$

## Installation
Run:
```bash
pip install numpy scipy matplotlib

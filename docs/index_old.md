---
layout: default
title: 1D Ising model simulations
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.js">
</script>

# 1D Ising Model with Transverse Field

The Hamiltonian of the system is:

$$H = - \sum_{n=1}^{N-1} \sigma^z_n \sigma^z_{n+1} - g \sum_{n=1}^{N} \sigma^x_n$$

where:
- $$N$$ is the number of spin-1/2 particles.
- $$g$$ is the transverse magnetic field strength.
- $$\sigma^x_n$$, $$\sigma^z_n$$ spin-1/2 operators associated with site $$n$$ in the direction of $$\hat{x}$$, $$\hat{z}$$ respectively

The model is in a **paramagnetic phase**, when $$g>>1$$. In this case ground state of the system has no degeneracy and it is a state in which all spins are pointing along $$\hat{x}$$:

$$\ket{GS}=\ket{\rightarrow\rightarrow\rightarrow...\rightarrow\rightarrow\rightarrow}$$

On the other hand, the most excited state in this phase would be the following one:

$$\ket{\leftarrow\leftarrow\leftarrow...\leftarrow\leftarrow\leftarrow}$$

The model is in a \textbf{ferromagnetic phase}, when $$g<<1$$ or $$g \to 0$$. A ground state of the system has a two-fold degeneracy and it can be one of the following states (or their superposition):

$$\ket{\uparrow\uparrow\uparrow...\uparrow\uparrow\uparrow}, \quad \ket{\downarrow\downarrow\downarrow...\downarrow\downarrow\downarrow}$$

On the other hand, the most excited state in this phase would be:

$$\ket{\uparrow\downarrow\uparrow...\downarrow\uparrow\downarrow} \quad \mathrm{or} \quad \ket{\downarrow\uparrow\downarrow...\uparrow\downarrow\uparrow}$$

We can make Jordan-Wigner transformation to spinless non-interacting fermions $$c_n$$, $$c^{\dagger}_n$$ given by:

$$\sigma^x_n = 1 - 2 c^{\dagger}_n c_n, \quad \sigma^z_n = (c_n + c^{\dagger}_n) \prod (1 - 2 c^{\dagger}_m c_m)$$

Our Hamiltonian becomes quadratic in creation and annihilation fermionic operators $$c_n$$, $$c^{\dagger}_n$$:

$$H = 2 \sum_{n=1}^{N} {g_n c^{\dagger}_n c_n} - 2 \sum (g_n c^{\dagger}_n c_n c_{n+1} c_n)$$ 

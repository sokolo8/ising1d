# 🧠 1D Transverse-Field Ising Model – Exact Simulation via Free Fermions

This repository contains a fast and accurate simulation of the **1D Ising model with a transverse field**, exploiting its exact solvability through the **Jordan-Wigner transformation** to free fermions.

---

## 🧩 Model Overview

The Hamiltonian of the system is:

$$
H = -\sum_{n=1}^{N-1} \sigma^z_n \sigma^z_{n+1} - \sum_{n=1}^{N} g_n \sigma^x_n
$$

- $N$: Number of spin-1/2 sites  
- $g_n$: Transverse magnetic field (can be time-dependent or spatially varying)  
- $\sigma^x_n, \sigma^z_n$: Pauli matrices acting on site $n$

---

## 🔁 Exact Solvability

Using the **Jordan-Wigner transformation**, the model is mapped to a system of **non-interacting fermions**:

$$
\sigma^x_n = 1 - 2 c^{\dagger}_n c_n,
$$
$$
\sigma^z_n = -(c_n + c^{\dagger}_n) \Pi_{mn} (1 - 2 c^{\dagger}_m c_m)
$$

This allows the Hamiltonian to be rewritten as a **quadratic fermionic form** that can be diagonalized via a **Bogoliubov transformation**.

---

## 🚀 Features

- ✅ Exact simulation of quantum quenches in 1D Ising model
- ✅ Time evolution using BdG equations
- ✅ Kibble-Zurek scaling analysis
- ✅ Suzuki-Trotter decomposition support
- ✅ Supports both static and time-dependent fields $g_n(t)$

---

## 📂 Contents

- `ising_model.py`: Core simulation engine
- `simulation.ipynb`: Interactive notebook to explore dynamics
- `docs/`: Optional LaTeX or Sphinx-based documentation
- `plots/`: Generated plots of excitation densities, etc.

---

## 📈 Sample Output

*(Insert plot or animation here if available)*

---

## 📄 Documentation

📓 **Notebook version**: [Open on GitHub](./README.ipynb)

---

## 📜 License

MIT License — free to use, cite, and modify. Feel free to build on it!

---

## 🤝 Acknowledgements

- Inspired by classic results on quantum phase transitions in 1D
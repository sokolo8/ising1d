## 1D Ising Model with Transverse Field
The Hamiltonian of the system is:

$$H = - \sum_{n=1}^{N-1} \sigma^z_n \sigma^z_{n+1} - \sum_{n=1}^{N} g_n\sigma^x_n$$

where:
- $N$ is the number of spin-1/2 particles.
- $g$ is the transverse magnetic field strength.
- $\sigma^x_n$, $\sigma^z_n$ are spin-1/2 operators at site $n$ in $\hat{x}, \hat{z}$ direction.

### Phases:
**Paramagnetic Phase** ($g \gg 1$):
- Ground state: $$\ket{GS}=\ket{\rightarrow\rightarrow\cdots\rightarrow}$$
- Most excited: $$\ket{\leftarrow\leftarrow\cdots\leftarrow}$$

**Ferromagnetic Phase** ($g \ll 1$):
- Degenerate ground states: $$\ket{\uparrow\uparrow\cdots\uparrow}, \quad \ket{\downarrow\downarrow\cdots\downarrow}$$
- Excited states: $$\ket{\uparrow\downarrow\uparrow\cdots}, \quad \ket{\downarrow\uparrow\downarrow\cdots}$$

We use the Jordan-Wigner transformation:
$$\sigma^x_n = 1 - 2 c^{\dagger}_n c_n, \quad \sigma^z_n = -(c_n + c^{\dagger}_n) \prod_{m<n} (1 - 2 c^{\dagger}_m c_m)$$

Resulting Hamiltonian:
$$H = 2 \sum_{n=1}^{N} g_n c^{\dagger}_n c_n - 2 \sum_{n=1}^{N-1} (g_n c^{\dagger}_n c_n + c_{n+1} c_n)$$

Diagonalized via Bogoliubov transformation:
$$H = \sum_{m=1}^{N}\omega_m \left(\gamma_m^{\dagger} \gamma_m - \frac{1}{2} \right)$$
$$c_n = \sum_{m=1}^{N} ( u_{n,m} \gamma_m + v_{n,m}^* \gamma^{\dagger}_m )$$

Stationary BdG equations:
$$\omega_m u^{\pm}_{n,m} = 2 g_n u^{\mp}_{n,m} - 2 u^{\mp}_{n\mp1, m}$$

Time-dependent BdG equation:
$$i \hbar \partial_{t} u^{\pm}_{n,m} = 2 g_n(t) u^{\mp}_{n,m} - 2 u^{\mp}_{n\mp1, m}$$
Boundary: $$u^{\pm}_{n+1}=u^{\pm}_{0}=0$$

Vector form:
$$i \hbar \partial_t \begin{pmatrix} \vec{u}^{+} \\ \vec{u}^{-} \end{pmatrix} = H_1(t)\begin{pmatrix} \vec{u}^{+} \\ \vec{u}^{-} \end{pmatrix}+H_2 \begin{pmatrix} \vec{u}^{+} \\ \vec{u}^{-} \end{pmatrix}$$

Formal solution:
$$\begin{pmatrix} \vec{u_t}^{+} \\ \vec{u_t}^{-} \end{pmatrix} = e^{-i \hbar (H_1(t) + H_2)}\begin{pmatrix} \vec{u_{t=0}}^{+} \\ \vec{u_{t=0}}^{-} \end{pmatrix}$$
Trotterized evolution:
$$U(t) \approx \prod_{n=0}^{N-1} e^{-i \hbar dt (H_1(ndt+dt/2) + H_2)}$$
$$e^{-i \hbar dt (H_1 + H_2)} = e^{-i \hbar \frac{dt}{2} H_2} e^{-i \hbar dt H_1} e^{-i \hbar \frac{dt}{2} H_2} + O(dt^3)$$

Hamiltonians in diagonal bases:
$$i \hbar \partial_t u_{n,m} = 2g_n(t) u_{n,m} + \dots$$
$$i \hbar \partial_t v_{n,m} = -2g_n(t) v_{n,m} + \dots$$
$$i \hbar \partial_t \tilde{u}^+_{n,m} = \dots - 2\tilde{u}^+_{n,m}$$
$$i \hbar \partial_t \tilde{u}^-_{n,m} = \dots + 2\tilde{u}^-_{n,m}$$

Where:
$$\tilde{u}^+_n = \frac{u^{+}_n + u^-_{n-1}}{\sqrt{2}}, \quad \tilde{u}^-_n = \frac{u^{+}_{n+1} - u^-_{n}}{\sqrt{2}}$$

Excitation density:
$$p = \frac{1}{N}\sum_{m}\bra{\psi(t_f)}\gamma^{\dagger}_m\gamma_m \ket{\psi(t_f)}$$

Where:
$$\gamma_m = \sum_{n=1}^{N} ( u^*_{n,m} c_n + v^*_{n,m} c^{\dagger}_n )$$
$$c_n = \sum_{m=1}^{N} ( u_{n,m}(t_f) \tilde{\gamma}_m + v_{n,m}^*(t_f) \tilde{\gamma}^{\dagger}_m )$$
with $\tilde{\gamma}_m \ket{\psi(t_f)}=0$

## Kibble-Zurek Mechanism
Spontaneous symmetry breaking creates topological defects:
$$\mathbb{Z}_2 \rightarrow \mathbb{1}$$

Near the critical point:
$$g \sim \left(\frac{t}{\tau_Q}\right)^r$$

Defect density:
$$d \sim \tau_Q^{-\frac{r\nu}{rz\nu+1}}$$
For Ising model: $z=\nu=1$, $r=2 \Rightarrow d \sim \tau_Q^{-2/3}$

### Example 1
Sine ramp from $g=2$ to $g_c=1$:
$$g(t) =
\begin{cases}
2, & \text{if } \frac{t}{\tau_Q} < -\frac{\pi}{2} \\
3/2 - 1/2 \sin(\frac{t}{\tau_Q}), & \text{if } -\frac{\pi}{2} \leq \frac{t}{\tau_Q} \leq \frac{\pi}{2} \\
1, & \text{if } \frac{t}{\tau_Q} > \frac{\pi}{2}
\end{cases}$$

Data for $L=100, 200, 300, 400, 500, 1000$ confirm power-law scaling.

### Example 2
Dispersion for infinite 1D Ising model:
$$\omega_k=2\sqrt{(g-\cos{k})^2+\sin^2{k}}$$
Near $k=0$:
$$\omega_k \approx 2k \Rightarrow c=\frac{\mathrm{d} \omega}{\mathrm{d}k}=2$$

Interpretation: if ramp speed is slower than $c$, parts of the system can communicate orientation across the transition, reducing excitations.

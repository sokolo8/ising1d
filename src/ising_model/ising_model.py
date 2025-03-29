import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh
from tqdm import tqdm # progressbar
from scipy.sparse import SparseEfficiencyWarning

import warnings
warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)

class Ising_Error(Exception):

    def __init__(self, message=None):
        super().__init__(message)

class Ising:

    def __init__(self, N, g_i, g_f, boundary='OBC'):
        self.N = N
        self.g_i = g_i
        self.g_f = g_f
        self.boundary = boundary
        self._initialize_hamiltonian()
        self._compute_initial_state()
        self._compute_final_state()

    def _initialize_hamiltonian(self):

        N = self.N
        self.External_field_Matrix = csc_matrix((2 * N, 2 * N))
        self.Basis_EFM = csc_matrix((2 * N, 2 * N))
        self.Couplings_Matrix_X = csc_matrix((2 * N, 2 * N))
        self.Basis_CMX = csc_matrix((2 * N, 2 * N))

        for i in range(2 * N):
            self.External_field_Matrix[i, 2 * N - i - 1] = 2.0

        for i in range(N):
            self.Basis_EFM[i, i] = 1 / np.sqrt(2)
            self.Basis_EFM[i + N, i + N] = -1 / np.sqrt(2)
            self.Basis_EFM[i, 2 * N - i - 1] = 1 / np.sqrt(2)
            self.Basis_EFM[i + N, N - i - 1] = 1 / np.sqrt(2)

        if self.boundary=='OBC':

            for i in range(0, N-1):
                self.Couplings_Matrix_X[i, 2 * N - i - 2] = -2.0
                self.Couplings_Matrix_X[N + i, N - i - 2] = -2.0

            for i in range(0, N-1):
                self.Basis_CMX[i, 2 * N - 2 * i - 3] = 1 / np.sqrt(2)
                self.Basis_CMX[i, 2 * N - 2 * i - 2] = 1 / np.sqrt(2)
                self.Basis_CMX[N + i, 2 * i + 1] = 1 / np.sqrt(2)
                self.Basis_CMX[N + i, 2 * i + 2] = -1 / np.sqrt(2)

            self.Basis_CMX[N - 1, 0] = 1
            self.Basis_CMX[2 * N - 1, 2 * N - 1] = 1

        if self.boundary=='PBC':

            for i in range(0, N-1):
                self.Couplings_Matrix_X[i, 2 * N - i - 2] = -2.0
                self.Couplings_Matrix_X[N + i, N - i - 2] = -2.0

            self.Couplings_Matrix_X[N - 1, 2 * N - 1] = -2.0
            self.Couplings_Matrix_X[2 * N - 1, N - 1] = -2.0

            for i in range(0, N-1):
                self.Basis_CMX[i, 2 * N - 2 * i - 2] = 1 / np.sqrt(2)
                self.Basis_CMX[i, 2 * N - 2 * i - 1] = 1 / np.sqrt(2)
                self.Basis_CMX[N + i, 2 * i + 2] = 1 / np.sqrt(2)
                self.Basis_CMX[N + i, 2 * i + 3] = -1 / np.sqrt(2)

            self.Basis_CMX[N - 1, 0] = 1 / np.sqrt(2)
            self.Basis_CMX[N - 1, 1] = 1 / np.sqrt(2)
            self.Basis_CMX[2 * N - 1, 0] = 1 / np.sqrt(2)
            self.Basis_CMX[2 * N - 1, 1] = -1 / np.sqrt(2)

        # Compute eigenvalues of Couplings matrix

        Diagonal_CMX = self.Basis_CMX.T @ self.Couplings_Matrix_X @ self.Basis_CMX
        self.Eigenvalues_of_CMX = Diagonal_CMX.diagonal()

        # Compute eigenvalues of External field matrix

        Diagonal_EFM = self.Basis_EFM.T @ self.External_field_Matrix @ self.Basis_EFM
        self.Eigenvalues_of_EFM = Diagonal_EFM.diagonal()

    def _compute_initial_state(self):

        # Compute the initial state in {u^-, u^+} representation

        self.Initial_H = self.External_field_Matrix * self.g_i + self.Couplings_Matrix_X
        self.initial_spectrum, self.initial_state = eigsh(self.Initial_H, k=self.N, which='LA')

    def _compute_final_state(self):

        # Compute the final state in {u^-, u^+} representation

        self.final_H = self.External_field_Matrix * self.g_f + self.Couplings_Matrix_X
        self.final_spectrum, self.final_state = eigsh(self.final_H, k=self.N, which='LA')

        # Compute the final state in {u, v} representation

        final_state_uv = self.Basis_EFM @ self.final_state

        self.U_final = final_state_uv[self.N - 1::-1, :]
        self.V_final = final_state_uv[self.N:2 * self.N, :]

    def homo_evolution(self, TTime, dt):

        Pi = np.pi
        Num = round(TTime / dt)

        g_i = self.g_i
        g_f = self.g_f
        N = self.N

        while Num < 1000:
            dt = dt / 2
            Num = round(TTime / dt)

        def g(t):
            if -Pi / 2.0 <= Pi / 2.0 - dt * t * Pi / TTime <= Pi / 2.0:
                return 0.5 * (g_i + g_f) + 0.5 * (g_i - g_f) * np.sin(Pi / 2.0 - dt * t * Pi / TTime)
            elif Pi / 2.0 - dt * t * Pi / TTime >= Pi / 2.0:
                return g_i
            elif Pi / 2.0 - dt * t * Pi / TTime <= -Pi / 2.0:
                return g_f


        Evolution_of_CMX = csc_matrix(np.diag(np.exp(-0.5j * dt * self.Eigenvalues_of_CMX)))

        # Construct the evolution matrices as sparse matrices
        evolution_matrix_right = self.Basis_EFM.T @ self.Basis_CMX @ Evolution_of_CMX
        evolution_matrix_left = Evolution_of_CMX @ self.Basis_CMX.T @ self.Basis_EFM

        # prepare the state in the convenient basis
        evolved_state = self.Basis_CMX.T @ self.initial_state

        for step in tqdm(range(Num)):
            # Calculate the complex phase
            magnetic_field_values = [g(step + 1) for _ in range(N)]
            rotation_angles = -2j * dt * np.flip(magnetic_field_values)
            concatenated_rotation_angles = np.concatenate((rotation_angles, -np.flip(rotation_angles)))
            exponential_rotation_angles = np.exp(concatenated_rotation_angles)

            # Apply the complex phase as a diagonal matrix
            diagonal_matrix_rotation_angles = csc_matrix(np.diag(exponential_rotation_angles))
            evolved_state = evolution_matrix_left @ diagonal_matrix_rotation_angles @ evolution_matrix_right @ evolved_state

        # Switch to {u^-, u^+} representation
        evolved_state_uv = self.Basis_EFM.T @ self.Basis_CMX @ evolved_state

        U_evolved = evolved_state_uv[N - 1::-1, :]
        V_evolved = evolved_state_uv[N:2 * N, :]

        term1 = np.transpose(U_evolved) @ self.V_final @ np.transpose(np.conjugate(self.V_final)) @ np.conjugate(U_evolved)
        term2 = np.transpose(U_evolved) @ self.V_final @ np.transpose(np.conjugate(self.U_final)) @ np.conjugate(V_evolved)
        term3 = np.transpose(V_evolved) @ self.U_final @ np.transpose(np.conjugate(self.V_final)) @ np.conjugate(U_evolved)
        term4 = np.transpose(V_evolved) @ self.U_final @ np.transpose(np.conjugate(self.U_final)) @ np.conjugate(V_evolved)

        sum_matrix = term1 + term2 + term3 + term4

        Number_of_kinks = sum_matrix.diagonal().sum()
        Density_of_kinks = Number_of_kinks / N

        Energy = sum(sum_matrix.diagonal() * self.final_spectrum) - 1/2 * sum(self.final_spectrum)
        Energy_GS = - 1/2 * sum(self.final_spectrum)
        Delta_Energy_per_bond = (Energy - Energy_GS) / N

        return np.real(Density_of_kinks), np.real(Delta_Energy_per_bond)
    
    def inhomo_evolution(self, TTime, dt, alpha):

        if self.boundary!='OBC':
            raise Ising_Error('You can only use it for an open chain')

        g_i = self.g_i
        g_f = self.g_f
        N = self.N

        Pi = np.pi

        Num = round(TTime / dt)

        while Num < 1000:
            dt = dt / 2
            Num = round(TTime / dt)

        R = (N - 1) / 2
        coords = [i - R for i in range(N)]
        v = Pi / (TTime * alpha) + np.sqrt(2) * R / TTime

        def g(t, n):
            if -Pi / 2.0 <= Pi / 2.0 + (np.abs(coords[n - 1]) - dt * v * t) * alpha <= Pi / 2.0:
                return 1 / 2 * (g_i + g_f) + 1 / 2 * (g_i - g_f) * np.sin(
                    Pi / 2.0 + (np.abs(coords[n - 1]) - dt * v * t) * alpha)
            elif Pi / 2.0 + (np.abs(coords[n - 1]) - dt * v * t) * alpha >= Pi / 2.0:
                return g_i
            elif Pi / 2.0 + (np.abs(coords[n - 1]) - dt * v * t) * alpha <= -Pi / 2.0:
                return g_f


        Evolution_of_CMX = csc_matrix(np.diag(np.exp(-0.5j * dt * self.Eigenvalues_of_CMX)))

        # Construct the evolution matrices as sparse matrices
        evolution_matrix_right = self.Basis_EFM.T @ self.Basis_CMX @ Evolution_of_CMX
        evolution_matrix_left = Evolution_of_CMX @ self.Basis_CMX.T @ self.Basis_EFM

        # prepare the state in the convenient basis
        evolved_state = self.Basis_CMX.T @ self.initial_state

        for step in tqdm(range(Num)):
            # Calculate the complex phase
            magnetic_field_values = [g(step + 1, j + 1) for j in range(N)]
            rotation_angles = -2j * dt * np.flip(magnetic_field_values)
            concatenated_rotation_angles = np.concatenate((rotation_angles, -np.flip(rotation_angles)))
            exponential_rotation_angles = np.exp(concatenated_rotation_angles)

            # Apply the complex phase as a diagonal matrix
            diagonal_matrix_rotation_angles = csc_matrix(np.diag(exponential_rotation_angles))
            evolved_state = evolution_matrix_left @ diagonal_matrix_rotation_angles @ evolution_matrix_right @ evolved_state

        # Switch to {u^-, u^+} representation
        evolved_state_uv = self.Basis_EFM.T @ self.Basis_CMX @ evolved_state

        U_evolved = evolved_state_uv[N - 1::-1, :]
        V_evolved = evolved_state_uv[N:2 * N, :]

        term1 = np.transpose(U_evolved) @ self.V_final @ np.transpose(np.conjugate(self.V_final)) @ np.conjugate(U_evolved)
        term2 = np.transpose(U_evolved) @ self.V_final @ np.transpose(np.conjugate(self.U_final)) @ np.conjugate(V_evolved)
        term3 = np.transpose(V_evolved) @ self.U_final @ np.transpose(np.conjugate(self.V_final)) @ np.conjugate(U_evolved)
        term4 = np.transpose(V_evolved) @ self.U_final @ np.transpose(np.conjugate(self.U_final)) @ np.conjugate(V_evolved)

        sum_matrix = term1 + term2 + term3 + term4

        Number_of_kinks = sum_matrix.diagonal().sum()
        Density_of_kinks = Number_of_kinks / N

        Energy = sum(sum_matrix.diagonal() * self.final_spectrum) - 1/2 * sum(self.final_spectrum)
        Energy_GS = - 1/2 * sum(self.final_spectrum)
        Delta_Energy_per_bond = (Energy - Energy_GS) / N

        return np.real(Density_of_kinks), np.real(Delta_Energy_per_bond)
    

    def compute_energy_per_bond_GS(self):
        """Computes observables per bond of final Hamiltonian at g = g_f"""

        N = self.N
        g_f = self.g_f

        U = self.U_final
        V = self.V_final

        term1 = V @ V.H
        term2 = U @ V.H

        E_ZZ = []

        E_ZZ.append(- 2 * term1[0, 1] - 2 * term2[1, 0])

        for i in range(1, N - 2):
            E_ZZ.append(- 2 * term1[i, i+1] - 2 * term2[i+1, i])

        E_ZZ.append(- 2 * term1[N-2, N-1] - 2 * term2[N-1, N-2])

        E_X = []

        for i in range(N):
            E_X.append(2 * g_f * term1[i, i] - g_f)

        E = []

        E.append(E_X[0] + E_X[1] / 2 + E_ZZ[0])

        for i in range(1, N - 2):
            E.append(E_X[i] / 2 + E_X[i+1] / 2 + E_ZZ[i])

        E.append(E_X[N-2] / 2 + E_X[N-1] + E_ZZ[N-2])

        return E, E_ZZ, E_X
    
    def compute_energy_GS(self):
        """Computes energy of final Hamiltonian at g = g_f"""

        N = self.N
        Energy = -1/2 * sum(self.final_spectrum) / N

        return Energy

if __name__ == '__main__':

    g_i = 2
    g_f = 1

    N = 100

    TTime = 100
    dt = 0.01
    alpha = 1/8

    system = Ising(N, g_i, g_f)

    print(system.compute_energy_GS())

    print(system.homo_evolution(TTime, dt))
    print(system.inhomo_evolution(TTime, dt, alpha))
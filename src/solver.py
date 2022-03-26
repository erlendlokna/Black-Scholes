import numpy as np
 
#function for creating a tridiag-coef matrix:
def tridiag(diags):
    M = len(diags[1])
    A = np.zeros((M, M))
    ks = [-1, 0, 1]
    for i in range(len(diags)):
        diag, k = diags[i], ks[i]
        A += np.diag(diag, k)
    return A

#abstract solver. Contains common features among the different methods
class Solver:
    def __init__(self, k_factor=1):  # Default empty computational grid
        self.k_factor = k_factor
        self.ibvp = None  # Initial/boundary value problem

        self.M = 1
        self.h = 0
        self.m = np.empty(self.M+1)
        self.x = np.empty(self.M+1)

        self.N = 1
        self.k = 0
        self.n = np.empty(self.N+1)
        self.t = np.empty(self.N+1)

    # Set computational grid {
    def set_grid(self, ibvp, M):
        self.ibvp = ibvp
        
        self.M = M
        self.h = ibvp.R / self.M
        self.m = np.arange(0, self.M+1)
        self.x = self.h * self.m
        
        self.N = int(self.M_to_N())
        self.k = self.k_factor * self.ibvp.T / self.N
        #self.k = (self.k_factor / 2) * (self.h/(ibvp.sigma*ibvp.R))**2
        self.n = np.arange(0, self.N+1)
        self.t = self.k * self.n
        
        #print("N =", self.N)
    # }

    def M_to_N(self):
        return int(np.ceil((2 * self.M**2 * self.ibvp.R**2 * self.ibvp.T * self.ibvp.sigma**2)))

    # Compute coefficients {
    def theta(self, m):
        return (self.k * (m * self.ibvp.sigma)**2)/2
        
    def phi(self, m):
        return (self.k * self.ibvp.r * m)/2
    
    def gamma(self, m):
        return 2 * self.theta(m) - self.k * self.ibvp.c
    
    def get_coeffs(self):
        raise NotImplementedError()
    # }

    # Set coefficient, addition and initial solution matrices {
    def get_A(self, coeffs):  # (M-1) x (N+1) Coefficient matrix
        return tridiag([coeffs[0][2:-1], coeffs[1][1:-1], coeffs[2][1:-2]])

    def get_B(self, coeffs):  # (M-1) x (N+1) Addition matrix
        B = np.zeros((self.M - 1, self.N + 1))
        top, bottom = self.ibvp.boundary(self.t)
        B[0, :] = coeffs[0][1] * top
        B[-1, :] = coeffs[2][-2] * bottom
        return B

    def get_U(self):  # (M+1) x (N+1) Solution matrix
        U = np.empty((self.M + 1, self.N + 1))
        U[:, 0] = self.ibvp.initial(self.x)  # initial condition at t=0
        top, bottom = self.ibvp.boundary(self.t)
        U[0, :] = top  # boundary condition at x=0
        U[-1, :] = bottom  # boundary condition at x=R
        return U
    # }

    # Compute solution matrix {
    def solve(self, ibvp, M):
        raise NotImplementedError()
    # }


class ForwardEuler(Solver):
    def __init__(self, k_factor=1):
        super().__init__(k_factor)
        self.name = "FE"
        
    def get_coeffs(self, m):
        theta, phi, gamma = self.theta(m), self.phi(m), self.gamma(m)
        return [theta - phi, 1 - gamma, theta + phi]

    def step(self, A, U, B, n):  # U^1 = A * U^0 + B
        return A @ U[1:-1, n] + B[:, n]
    
    def solve(self, ibvp, M):
        self.set_grid(ibvp, M)

        coeffs = self.get_coeffs(self.m)
        A = self.get_A(coeffs)
        B = self.get_B(coeffs)
        U = self.get_U()
        for n in range(0, self.N):
            U[1:-1, n+1] = self.step(A, U, B, n)
        return self.x, self.t, U

    
class BackwardEuler(Solver):
    def __init__(self, k_factor=1):
        super().__init__(k_factor)
        self.name = "BE"

    def get_coeffs(self, m):
        theta, phi, gamma = self.theta(m), self.phi(m), self.gamma(m)
        return [-(theta - phi), 1 + gamma, -(theta + phi)]

    def step(self, A, U, B, n):  # U^1 = A * U^0 + B
        return np.linalg.solve(A, U[1:-1, n] + B[:, n+1])
    
    def solve(self, ibvp, M):
        self.set_grid(ibvp, M)

        # Matrix equation in each time step: A * U^1 = U^0 + B
        coeffs = self.get_coeffs(self.m)
        A = self.get_A(coeffs)
        B = self.get_B(coeffs)
        U = self.get_U()
        for n in range(0, self.N):
            U[1:-1, n+1] = self.step(A, U, B, n)
        return self.x, self.t, U
    


class CrankNicolson(Solver):
    def __init__(self, k_factor=1):
        super().__init__(k_factor)
        self.name = "CN"

    def get_coeffs(self, m):
        theta, phi, gamma = self.theta(m), self.phi(m), self.gamma(m)
        return [
            [-(theta - phi)/2, 1 + gamma/2, -(theta + phi)/2],
            [(theta - phi)/2, 1 - gamma/2, (theta + phi)/2]
        ]

    def step(self, A1, A0, U, B, n):
        return np.linalg.solve(A1, A0 @ U[1:-1, n] + B)

    def solve(self, ibvp, M):
        self.set_grid(ibvp, M)

        # Matrix equation in each time step: A1 * U^1 = A0 * U^0 + B
        A1_coeffs, A0_coeffs = self.get_coeffs(self.m)
        A1, A0 = self.get_A(A1_coeffs), self.get_A(A0_coeffs)
        B1, B0 = self.get_B(A1_coeffs), self.get_B(A0_coeffs)
        U = self.get_U()
        for n in range(0, self.N):
            B = B1[:, n+1] + B0[:, n]
            U[1:-1, n+1] = self.step(A1, A0, U, B, n)
        return self.x, self.t, U
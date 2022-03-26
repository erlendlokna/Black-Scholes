# IBVP = initial/boundary value problem
# NotImplementedError() = abstract method

######  Initial and boundary conditions ###### 

import numpy as np


class IBVP:
    def __init__(self, 
                 sigma, r, c, 
                 R, K, T):
        self.sigma, self.r, self.c = sigma, r, c
        self.R, self.K, self.T = R, K, T
    
    def boundary(self, t):
        left = np.zeros_like(t)
        right = self.initial(self.R) * np.ones_like(t)
        return left, right
        
    def initial(self, x):
        raise NotImplementedError()
        
        
class EUput(IBVP):
    def __init__(self, 
                 sigma, r, c, 
                 R, K, T):
        super().__init__(sigma, r, c, R, K, T)
        self.name = "EUput"
    
    def boundary(self, t):
        return self.K * np.exp(-self.c * t), np.zeros_like(t)
            
    def initial(self, x):
        return np.maximum(self.K - x, 0)
            
            
class ButterflySpread(IBVP):
    def __init__(self, 
                 sigma, r, c, 
                 R, K, H, T):
        super().__init__(sigma, r, c, R, K, T)
        self.H = H
        self.name = "Butterfly Spread"
            
    def initial(self, x):
        K, H, m = self.K, self.H, np.maximum
        return m(x - K, 0) - 2*m(x - (K + H), 0) + m(x - (K + 2 * H), 0)
            

class BinaryCall(IBVP):
    def __init__(self, 
                 sigma, r, c, 
                 R, K, T):
        super().__init__(sigma, r, c, R, K, T)
        self.name = "Binary Call"
            
    def initial(self, x):
        return np.maximum(np.sign(x - self.K), 0)



class Exact(IBVP):
    #Object to hold the exact solution with its boundary conditions and initial values.
    def __init__(self, 
                 sigma, r, c, 
                 R, K, T):
        super().__init__(sigma, r, c, R, K, T)

    def boundary(self, t):
        return self.solution(0, t), np.zeros_like(t)
            
    def initial(self, x):
        return self.solution(x, 0)

    def solution(self, x, t):
        a = - np.pi**2 * np.sin(np.pi * x)
        b = 0.5 * self.sigma**2 * x**2 * np.pi**2 * np.sin(np.pi * x)
        c = self.r * x * np.pi * np.cos(np.pi * x)
        d = self.c * np.sin(np.pi * x)
        return a + b + c + d

#it is now a simple task to add the new IBVP based on the "exact" solution to the solver and calculate the error by comparing numerical and exact solution.
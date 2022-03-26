from src import boundary as boundary
from src import convergenceError as error
from src import graphic as graphic
from src import solver as solver

import matplotlib.pyplot as plt
import numpy as np

def solve(ibvp, SolverClass, M, k_factor=1, title=""):
    solver = SolverClass(k_factor)
    graphic.plot_solution(*solver.solve(ibvp, M), title + " " + solver.name + " " + ibvp.name)


def test_stability(ibvp, SolverClass, M, k_factor=3):
    solve(ibvp, SolverClass, M, 1, "stable")
    solve(ibvp, SolverClass, M, k_factor, "unstable")


def main():
    sigma, r, c = 1, .0, .1
    R, K, H, T = 1, .1, .2, 1
    M = 20

    eu_put = boundary.EUput(sigma, r, c, R, K, T)
    binary_call = boundary.BinaryCall(sigma, r, c, R, K, T)
    butterfly_spread = boundary.ButterflySpread(sigma, r, c, R, K, H, T)

    test_stability(eu_put, solver.ForwardEuler, M, 3)
    test_stability(binary_call, solver.ForwardEuler, M, 3)
    test_stability(butterfly_spread, solver.ForwardEuler, M, 3)


    #list containing different M̈́'s
    Ms = np.array([10, 15, 20, 25, 30, 35])

    #the exact solution
    _exact = boundary.Exact(sigma, r, c, R, K, T) 

    #calculating errors and visualising errors in a log log plot:
    norm = 2
    error.convergenceError(_exact, solver.ForwardEuler, Ms, norm = norm, txt="FE error") 
    error.convergenceError(_exact, solver.BackwardEuler, Ms,norm = norm, txt="BE error")
    error.convergenceError(_exact, solver.CrankNicolson, Ms, norm = norm, txt="CN error")

if __name__ == "__main__":
    main()


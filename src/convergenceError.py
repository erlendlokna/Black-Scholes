import numpy as np
import matplotlib.pyplot as plt

def convergenceError(ExactSolution, SolverClass, Ms, norm=2, k_factor=1, txt = "log log plot"):
    """
    Input: 
        ExactSolution: exact solution object with its boundary, initial values and solution functions.
        SolverClass: object to solve the PDE.
        Ms: list containing M's to be calculated for.
        norm: norm of the error.
    """
    solver = SolverClass() #initialise the object
    errors = np.zeros(len(Ms)) #errors during simulations
    Ns = np.zeros(len(Ms)) #N's list.

    for i in range(len(Ms)):        
        solver.set_grid(ExactSolution, Ms[i]) #updating grid
        
        N = solver.M_to_N() #calculating the N that satisfies the stability condition.

        current_errors = [] #list with current errors

        x, t, U = solver.solve(ExactSolution,  Ms[i]) #solving the PDE numerically

        for j in range(N):
            exact = ExactSolution.solution(t[j], x) #calculating exact solution
            current_errors.append(solver.h**(1/norm) * np.linalg.norm(np.abs(U[:, j] - exact), norm)) #calculating error

            range_x = np.linspace(0, x[-1], 200) #x-values to be interpolat

            interpolated_values = np.interp(range_x, x, U[:, j]) #piecevise linear interpolation. 
            exact_values = ExactSolution.solution(t[j], range_x)

            current_errors.append(max(np.array(interpolated_values) - np.array(exact_values)))

        Ns[i] = N #storing results.
        errors[i] = solver.k**(1/norm) * np.linalg.norm(current_errors, norm)

    fig = plt.figure(figsize=(6,5)) #visualising result.
    ax = fig.gca()
    ax.plot(np.log(Ns), np.log(errors), '*-', c='royalblue')
    ax.set_title(txt)
    ax.set_xlabel("log(N)")
    ax.set_ylabel("log(|error|)")
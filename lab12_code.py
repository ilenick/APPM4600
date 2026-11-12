
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time


def driver():

    '''Compare direct solve and LU decomposition methods
       for solving square linear systems Ax=b'''

    # Matrix sizes to test
    N_values = [100, 500, 1000, 2000, 4000, 5000]

    # Store timing results
    direct_times = []
    lu_factor_times = []
    lu_solve_times = []
    lu_total_times = []

    # Loop through each matrix size
    for N in N_values:

        print(f"\nTesting N = {N}")

        # Generate random matrix and right-hand side
        A = np.random.rand(N, N)
        b = np.random.rand(N, 1)

        # -----------------------------------------
        # Method 1: Direct solution using solve()
        # -----------------------------------------

        start = time.perf_counter()

        x = scila.solve(A, b)

        end = time.perf_counter()

        direct_time = end - start
        direct_times.append(direct_time)

        # Calculate residual
        r = la.norm(A @ x - b)

        print(f"Regular solve time: {direct_time:.6f} seconds")
        print(f"Regular solve residual: {r:.6e}")

        # -----------------------------------------
        # Method 2: LU decomposition
        # -----------------------------------------

        # (i) Construct LU factorization
        start = time.perf_counter()

        lu, piv = scila.lu_factor(A)

        end = time.perf_counter()

        factor_time = end - start
        lu_factor_times.append(factor_time)

        # (ii) Solve using LU factors
        start = time.perf_counter()

        x_lu = scila.lu_solve((lu, piv), b)

        end = time.perf_counter()

        solve_time = end - start
        lu_solve_times.append(solve_time)

        # Calculate total LU time
        total_time = factor_time + solve_time
        lu_total_times.append(total_time)

        # Calculate residual
        r_lu = la.norm(A @ x_lu - b)

        print(f"LU factorization time: {factor_time:.6f} seconds")
        print(f"LU solve time: {solve_time:.6f} seconds")
        print(f"LU total time: {total_time:.6f} seconds")
        print(f"LU solve residual: {r_lu:.6e}")

    # -----------------------------------------
    # Print summary table
    # -----------------------------------------

    print("\n\nTIMING SUMMARY")
    print("-" * 85)
    print(f"{'N':<10}{'Direct Solve':<20}{'LU Factor':<20}"
          f"{'LU Solve':<20}{'LU Total':<15}")
    print("-" * 85)

    for i, N in enumerate(N_values):
        print(f"{N:<10}{direct_times[i]:<20.6f}"
              f"{lu_factor_times[i]:<20.6f}"
              f"{lu_solve_times[i]:<20.6f}"
              f"{lu_total_times[i]:<15.6f}")

    # -----------------------------------------
    # Plot timing results
    # -----------------------------------------

    plt.figure(figsize=(10, 6))

    plt.loglog(N_values, direct_times, 'o-',
               label='Regular solve')

    plt.loglog(N_values, lu_factor_times, 's-',
               label='LU factorization')

    plt.loglog(N_values, lu_solve_times, '^-',
               label='LU solve')

    plt.loglog(N_values, lu_total_times, 'd--',
               label='LU total')

    plt.xlabel('Matrix size N')
    plt.ylabel('Time (seconds)')
    plt.title('Comparison of Linear System Solution Methods')

    plt.legend()
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    driver()
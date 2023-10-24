from pysat.formula import CNF, WCNF
from pysat.examples.fm import FM
from ode import ODESolver
import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal, getcontext
getcontext().prec = 1

class Solver:
    def __init__(self, file_name: str, dt=0.01, T=10):
        self.dt = dt
        self.T = T

        f1 = CNF(from_file=file_name)
        self.c = np.zeros((len(f1.clauses), f1.nv))
        for i, clause in enumerate(f1.clauses):
            for var in clause:
                self.c[i][abs(var) - 1] = (var > 0) * 2 -1
        
        print(self.c)
        
        # M clauses, I variables
        self.M, self.I = self.c.shape

        print(f"M: {self.M}, I: {self.I}")

        # self.initial_a = np.random.normal(0,1, self.M)
        self.initial_a = np.ones(self.M) * 1
        self.initial_s = np.random.uniform(-1,1, self.I)
        self.initial = np.concatenate((self.initial_s, self.initial_a))

        self.solver = ODESolver(self.derivative, precision=1)

    def K(self, s, a, m, i=None):
        mask = np.ones(self.I, dtype=bool)
        if i is not None:
            mask[i] = False
        return np.prod(1 - self.c[m][mask] * s[mask])
    
    def ds_i(self, s, a, i):
        return sum([2*a[m]*self.c[m][i]*self.K(s, a, m, i) * self.K(s, a, m) for m in range(self.M)])

    def da_m(self, s, a, m):
        return a[m] * self.K(s, a, m)
    
    def V(self, s, a):
        K_values = self.K(s, a, np.arange(self.M))  # Calculate K values for all m
        return np.sum(a * K_values**2)

    def derivative(self, t, state):
        s, a = state[:self.I], state[self.I:]
        ds = np.array([self.ds_i(s, a, i) for i in range(len(s))])
        # da = np.array([self.da_m(s, a, m) for m in range(len(a))])
        da = np.ones(self.M)

        return np.concatenate((ds, da))
    
    def integrate(self):
        print("SOLVING....")
        times, values = self.solver.solve(self.initial, 0, 1, .01)
        return values[:,:self.I]

if __name__ == "__main__":
    soln = Solver('test_files/coloring.cnf')
    ts = soln.integrate()
    last = ts[-1]
    print(last)
    results = [index + 1 if value > 0 else -index - 1 for index, value in enumerate(ts[-1])]
    print(results)
    plt.ylim(-1.1, 1.1)
    plt.plot(ts)
    plt.legend(np.arange(soln.I)+1)
    plt.show()
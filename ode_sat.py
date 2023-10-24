import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

np.random.seed(0)

class OdeSat:

    def __init__(self, clauses: np.ndarray, resolution=1000, time=5*10**6):
        self.resolution = resolution
        self.time = time

        max_variable = max([max(np.abs(clause)) for clause in clauses])

        self.c = np.zeros((len(clauses), max_variable), dtype=int)
        for clause_index, clause in enumerate(clauses):
            for variable in clause:
                self.c[clause_index][abs(variable) - 1] = 1 if variable > 0 else -1

        print(self.c)
        
        # self.c = np.array([
        #     [1, 1, -1],
        #     [1, -1, 1],
        #     [1, -1, -1],
        #     [-1, 1, 1],
        #     [-1, 1, -1],
        #     [-1, -1, 1]
        # ])
        
        # M clauses, I variables
        self.M, self.I = self.c.shape

        print(f"M: {self.M}, I: {self.I}")

        self.initial_s = np.random.uniform(-1,1, self.I)
        self.initial_a = self.time * np.ones(self.M)
        self.initial = np.concatenate((self.initial_s, self.initial_a))

    def K(self, s, a, m, i = None):
        out = 1
        for j in range(self.I):
            if j == i: continue
            out *= 1-self.c[m][j]*s[j]
        # return 2**(-self.I + (i is not None))*out
        return out
    
    def ds_i(self, s, a, i):
        return (2**(-2 * self.I + 1))*sum([2*a[m]*self.c[m][i]*self.K(s, a, m, i) * self.K(s, a, m) for m in range(self.M)])

    def da_m(self, s, a, m):
        return a[m] * 2**(-self.I)*self.K(s, a, m)
    
    def V(self, s, a):
        return 2**(-self.I)*sum(a[m] * self.K(s, a, m)**2 for m in range(self.M))

    def derivative(self, state, t):
        s, a = state[:self.I], state[self.I:]
        ds = np.array([self.ds_i(s, a, i) for i in range(len(s))])
        da = np.array([self.da_m(s, a, m) for m in range(len(a))])

        return np.concatenate((ds, da))

    def integrate(self):
        t = np.linspace(0, self.time, self.resolution)
        timeseries = odeint(self.derivative, self.initial, t)
        return timeseries[:,:self.I]

    def run(self):
        timeseries = self.integrate()
        return timeseries


if __name__ == "__main__":
    from pysat.formula import CNF
    f1 = CNF(from_file='test_files/coloring.cnf')
    ode_sat = OdeSat(clauses=f1.clauses)
    ts = ode_sat.integrate()
    print(ts)
    last = ts[-1]
    print(last.tolist())
    print([index + 1 if value > 0 else -index - 1 for index, value in enumerate(last)])
    plt.plot(ts)
    plt.ylim(-1.1, 1.1)
    plt.legend(np.arange(ode_sat.I)+1)
    plt.show()
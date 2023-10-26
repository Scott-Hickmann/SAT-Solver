import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

np.random.seed(0)


class OdeSat:
    def __init__(self, clauses: np.ndarray, resolution=10000, time=2):
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

        self.initial_s = np.random.uniform(-1, 1, self.I)
        self.initial_a = self.time * np.ones(self.M)
        self.initial = np.concatenate((self.initial_s, self.initial_a))

    def K(self, s, a, m, i=None):
        c_vals = np.array(self.c[m])
        if i is not None:
            c_vals[i] = 0
        out = np.prod(1 - c_vals * s)
        return out

    def ds_i(self, s, a, i):
        c_vals = np.array(self.c)[:, i]
        K_vals_without_i = np.array([self.K(s, a, m) for m in range(self.M)])
        K_vals_with_i = np.array([self.K(s, a, m, i) for m in range(self.M)])
        return np.sum(2 * a * c_vals * K_vals_with_i * K_vals_without_i)

    def da_m(self, s, a, m):
        return a[m] * self.K(s, a, m)

    def derivative(self, state, t):
        print(f"t: {t}")
        s, a = state[: self.I], state[self.I :]
        ds = np.array([self.ds_i(s, a, i) for i in range(len(s))])
        da = np.array([self.da_m(s, a, m) for m in range(len(a))])

        d = np.concatenate((ds, da))
        max = np.amax(np.absolute(d))
        if max > 100:
            raise Exception(
                f"Stiff ODE. Step size is too large, leading to dx of {max}"
            )
        return d

    def integrate(self):
        try:
            t = np.linspace(0, self.time, self.resolution)
            # timeseries = odeint(self.derivative, self.initial, t)

            dt = t[1] - t[0]
            timeseries = np.empty((len(t), len(self.initial)))
            timeseries[0] = self.initial
            for i in range(1, len(t)):
                timeseries[i] = timeseries[i - 1] + dt * self.derivative(
                    timeseries[i - 1], t[i - 1]
                )

            return timeseries[:, : self.I]
        except Exception:
            raise

    def run(self):
        timeseries = self.integrate()
        return timeseries


def test_stiffness(clauses, time_limit=1):
    res = 10
    while res < 100000:
        try:
            ode_sat = OdeSat(clauses=f1.clauses, resolution=res, time=time_limit)
            ts = ode_sat.integrate()
        except:
            print("SKIP")
            res = int(res * 1.5)
            continue
        return res


if __name__ == "__main__":
    from pysat.formula import CNF

    f1 = CNF(from_file="test_files/coloring_basic.cnf")
    res = test_stiffness(f1.clauses)
    print("RES:", res)
    ode_sat = OdeSat(clauses=f1.clauses, resolution=109, time=1)
    ts = ode_sat.integrate()
    print(ts)
    last = ts[-1]
    print(last.tolist())
    print([index + 1 if value > 0 else -index - 1 for index, value in enumerate(last)])
    plt.plot(ts)
    plt.ylim(-1.1, 1.1)
    plt.legend(np.arange(ode_sat.I) + 1)
    plt.savefig("out/output_usa_euler_2s_10000res.png")

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

np.random.seed(0)


class OdeSat:
    def __init__(self, clauses: np.ndarray, resolution=1000, time=4):
        print("GO", resolution, time)
        self.resolution = resolution
        self.time = time

        max_variable = max([max(np.abs(clause)) for clause in clauses])

        self.c = np.zeros((len(clauses), max_variable), dtype=int)
        for clause_index, clause in enumerate(clauses):
            for variable in clause:
                self.c[clause_index][abs(variable) - 1] = 1 if variable > 0 else -1

        # print(self.c)

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
        c_vals = np.array(self.c)[m]
        if i is not None:
            c_vals[i] = 0
        out = np.prod(1 - c_vals * s)
        return out

    def V(self, s, a):
        K_vals = np.array([self.K(s, a, m) for m in range(self.M)])
        return a @ (K_vals**2)

    def ds_i(self, s, a, i):
        c_vals = np.array(self.c)[:, i]
        K_vals_without_i = np.array([self.K(s, a, m) for m in range(self.M)])
        K_vals_with_i = np.array([self.K(s, a, m, i) for m in range(self.M)])
        return np.sum(2 * a * c_vals * K_vals_with_i * K_vals_without_i)

    def da_m(self, s, a, m):
        return a[m] * self.K(s, a, m)

    def derivative(self, state):
        # print(f"t: {t}")
        s, a = state[: self.I], state[self.I :]
        ds = np.array([self.ds_i(s, a, i) for i in range(len(s))])
        da = np.array([self.da_m(s, a, m) for m in range(len(a))])

        d = np.concatenate((ds, da))
        return d

    def integrate(self):
        try:
            t = np.linspace(0, self.time, self.resolution)

            # ODE Int
            # timeseries = odeint(self.derivative, self.initial, t)

            # Euler
            dt = t[1] - t[0]
            timeseries = np.empty((len(t), len(self.initial)))
            timeseries[0] = self.initial
            # print("t shape", timeseries.shape)
            for i in range(1, len(t)):
                timeseries[i] = timeseries[i - 1] + dt * self.derivative(
                    timeseries[i - 1]
                )

                if np.isnan(timeseries).any() or np.any(
                    (timeseries < -1.5) | (timeseries > 1.5)
                ):
                    # print("stiff")
                    print(np.isnan(timeseries).any())
                    raise Exception(f"Stiff ODE. Step size is too large.")

            sTs = timeseries[:, : self.I]
            aTs = timeseries[:, self.I :]

            return sTs, aTs
        except Exception as e:
            print("exception occurred", e)
            raise

    def run(self):
        timeseries = self.integrate()
        return timeseries


def test_stiffness(clauses, time_limit=1):
    res = 10
    while res < 10000000:
        print(res)
        try:
            np.random.seed(0)
            ode = OdeSat(clauses=clauses, resolution=res, time=time_limit)
            ts = ode.integrate()
        except:
            res = int(res * 1.2)
            continue
        print("FOUND WORKING RES:", res)
        return res
    return -1


if __name__ == "__main__":
    from pysat.formula import CNF

    f1 = CNF(from_file="test_files/coloring_basic.cnf")
    result_name = "coloring_unoptimized"

    res = test_stiffness(f1.clauses)
    print("resolution", res)

    np.random.seed(0)
    ode_sat = OdeSat(clauses=f1.clauses, resolution=res, time=1)
    sTs, aTs = ode_sat.integrate()

    print(sTs)
    last = sTs[-1]
    print(last.tolist())
    print([index + 1 if value > 0 else -index - 1 for index, value in enumerate(last)])
    print("Final resolution:", res)
    plt.plot(sTs)
    plt.ylim(-1.1, 1.1)
    plt.legend(np.arange(ode_sat.I) + 1)
    plt.savefig(f"out/{result_name}.png")

    plt.clf()
    vs = [ode_sat.V(s, a) for s, a in zip(sTs, aTs)]
    plt.plot(vs)
    plt.savefig(f"out/{result_name}_v.png")

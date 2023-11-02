import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

np.random.seed(0)


class StiffException(Exception):
    pass


class OdeSat:
    def __init__(self, clauses: np.ndarray, resolution=1000, time=4):
        self.resolution = resolution
        self.time = time

        max_variable = max([max(np.abs(clause)) for clause in clauses])

        self.c = np.zeros((len(clauses), max_variable), dtype=int)
        for clause_index, clause in enumerate(clauses):
            for variable in clause:
                self.c[clause_index][abs(variable) - 1] = 1 if variable > 0 else -1

        # print(self.c)
        self.cT = self.c.T

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

    def K(self, s, m):
        return np.prod(1 - self.c[m] * s)

    # def Ki(self, K_vals, s, m, i):
    #     return K_vals[m] / (1 - self.c[m, i] * s[i])

    def V(self, s, a, K_vals=None):
        if np.all(K_vals) == None:
            K_vals = np.array([self.K(s, m) for m in range(self.M)])
        return a @ (K_vals**2)

    # def ds_i(self, s, a, i, K_vals):
    #     c_vals = self.c[:, i]
    #     denominator = 1 - (c_vals * s[i])
    #     Ki_vals = K_vals / denominator
    #     return 2 * a @ (c_vals * Ki_vals * K_vals)

    # def da_m(self, K_vals, a, m):
    #     return a[m] * K_vals[m]

    def da(self, a, K_vals):
        return a * K_vals

    def ds(self, s, a, K_vals):
        s_broadcasted = s[:, np.newaxis]
        denominator = 1 - (self.cT * s_broadcasted)
        Ki_vals = np.where(denominator != 0, K_vals / denominator, 0)
        intermediate_result = self.cT * Ki_vals * K_vals
        result = 2 * np.sum(a * intermediate_result, axis=1)
        return result

    def derivative(self, state):
        # print(f"t: {t}")
        s, a = state[: self.I], state[self.I :]
        K_vals = np.array([self.K(s, m) for m in range(self.M)])
        ds = self.ds(s, a, K_vals)
        da = self.da(a, K_vals)

        return np.concatenate((ds, da)), K_vals

    def integrate(self):
        t = np.linspace(0, self.time, self.resolution)

        # ODE Int
        # timeseries = odeint(self.derivative, self.initial, t)

        # Euler
        dt = t[1] - t[0]
        timeseries = np.empty((len(t), len(self.initial)))
        timeseries[0] = self.initial
        # print("t shape", timeseries.shape)
        for i in range(1, len(t)):
            ds_da, k_vals = self.derivative(timeseries[i - 1])
            timeseries[i] = timeseries[i - 1] + dt * ds_da

            sTs = timeseries[i, : self.I]

            if np.isnan(sTs).any() or np.any((sTs < -1.2) | (sTs > 1.2)):
                # print("stiff")
                print(np.isnan(sTs).any())
                print("i:", i)
                print(np.nanmax(sTs, 0))
                raise StiffException(f"Stiff ODE. Step size is too large.")

            sTs = timeseries[:i, : self.I]
            aTs = timeseries[:i, self.I :]

            vs = [
                self.V(s, a, k_vals)
                for s, a in zip(np.expand_dims(sTs[-1], 0), np.expand_dims(aTs[-1], 0))
            ]
            if np.min(vs) < 1e-3:
                print(f"ending early, {i / len(t) * 100}% of max scheduled time")
                return sTs, aTs

        return sTs, aTs

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
        except StiffException as e:
            print("NEXT")
            res = int(res * 1.2)
            continue
        print("FOUND WORKING RES:", res)
        return res
    return -1


if __name__ == "__main__":
    from pysat.formula import CNF

    f1 = CNF(from_file="test_files/factor8.cnf")
    result_name = "factor8"

    res = test_stiffness(f1.clauses)
    print("resolution", res)

    np.random.seed(0)
    ode_sat = OdeSat(clauses=f1.clauses, resolution=res, time=1)
    sTs, aTs = ode_sat.integrate()

    last = sTs[-1]
    print(last.tolist())
    print([index + 1 if value > 0 else -index - 1 for index, value in enumerate(last)])
    # print("Final resolution:", res)
    plt.plot(sTs)
    plt.ylim(-1.1, 1.1)
    plt.legend(np.arange(ode_sat.I) + 1)
    plt.savefig(f"out/{result_name}.png")

    plt.clf()
    vs = [ode_sat.V(s, a) for s, a in zip(sTs, aTs)]
    plt.plot(vs)
    plt.savefig(f"out/{result_name}_v.png")

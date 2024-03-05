import numpy as np
from scipy.integrate import odeint, solve_ivp
from matplotlib import pyplot as plt
from solver import verify

np.random.seed(0)

KTs = []
aTs = []

solution = [
    -1,
    2,
    -3,
    -4,
    5,
    -6,
    -7,
    8,
    -9,
    -10,
    11,
    -12,
    13,
    -14,
    15,
    -16,
    17,
    -18,
    -19,
    -20,
    -21,
    -22,
]


class OdeSat:
    def __init__(self, clauses: np.ndarray, resolution, time):
        self.resolution = resolution
        self.time = time

        max_variable = max([max(np.abs(clause)) for clause in clauses])

        self.c = np.zeros((len(clauses), max_variable), dtype=int)
        for clause_index, clause in enumerate(clauses):
            for variable in clause:
                self.c[clause_index][abs(variable) - 1] = 1 if variable > 0 else -1

        print(self.c)
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

        self.scale_a = 1
        # self.initial_s = np.random.uniform(-1, 1, self.I)
        # self.initial_s = np.ones(self.I)
        self.initial_s = np.array(
            [-0.5 if solution[i] > 0 else 0.5 for i in range(self.I)], dtype=float
        )
        self.initial_a = self.time * np.ones(self.M)
        # self.initial_a = np.ones(self.M)
        self.initial = np.concatenate((self.initial_s, self.initial_a))

    def K(self, s, m):
        return np.min(1 - self.c[m] * s)
        # return np.prod(1 - self.c[m] * s)

    def Ki(self, s, m, i, K_vals):
        variable_evaled = 1 - self.c[m] * s
        variable_evaled[self.c[m] == 0] = np.inf
        min_index = np.argmin(variable_evaled)
        # if m == 19 and i == 3:
        #     print(self.c[m], min_index)
        #     print(variable_evaled)
        return 1 * (min_index == i)

    def V(self, s, a):
        K_vals = np.array([self.K(s, m) for m in range(self.M)])
        return a @ (K_vals**2)

    def E(self, s):
        K_vals = np.array([self.K(s, m) for m in range(self.M)])
        return K_vals @ K_vals

    def ds_i(self, s, a, i, K_vals):
        c_vals = self.c[:, i]
        Ki_vals = np.array([self.Ki(s, m, i, K_vals) for m in range(self.M)])
        res = np.sum(2 * a * c_vals * Ki_vals * K_vals)
        if i == 1:
            print(i, res)
            print(c_vals)
            print(Ki_vals)
            print(K_vals)
        return res

    # def da_m(self, K_vals, a, m):
    #     return a[m] * K_vals[m]

    def da(self, a, K_vals):
        return a * K_vals

    def ds(self, s, a, K_vals, t):
        return np.array([self.ds_i(s, a, i, K_vals) for i in range(self.I)])
        # s_broadcasted = s[:, np.newaxis]
        # denominator = 1 - (self.cT * s_broadcasted)
        # Ki_vals = np.where(denominator != 0, K_vals / denominator, 0)
        # intermediate_result = self.cT * Ki_vals * K_vals
        # result = 2 * np.sum(a * intermediate_result, axis=1)
        # return result

    def derivative(self, t, state):
        print(f"t: {t}")
        s, a = state[: self.I], state[self.I :]
        K_vals = np.array([self.K(s, m) for m in range(self.M)])
        ds = self.ds(s, a, K_vals, t)
        da = self.da(a, K_vals)
        # Normalize as between 1 and 2
        KTs.append(K_vals)
        aTs.append(a)

        d = np.concatenate((ds, da))
        max = np.amax(np.absolute(d))
        # if max > 100:
        #     raise Exception(
        #         f"Stiff ODE. Step size is too large, leading to dx of {max}"
        #     )
        return d

    def integrate(self):
        t = np.linspace(0, self.time, self.resolution)

        # ODE Int
        # timeseries = odeint(self.derivative, self.initial, t)

        # Solve IVP
        res = solve_ivp(
            self.derivative,
            (0, self.time),
            self.initial,
            method="LSODA",
        )
        timeseries = res.y.T

        # Euler
        # dt = t[1] - t[0]
        # timeseries = np.empty((len(t), len(self.initial)))
        # timeseries[0] = self.initial
        # for i in range(1, len(t)):
        #     next_vals = timeseries[i - 1] + dt * self.derivative(
        #         t[i - 1], timeseries[i - 1]
        #     )
        #     next_vals[: self.I] = np.clip(next_vals[: self.I], -1, 1)
        #     timeseries[i] = next_vals

        sTs = timeseries[:, : self.I]
        aTs = timeseries[:, self.I :]

        return sTs, aTs

    def run(self):
        timeseries = self.integrate()
        return timeseries


if __name__ == "__main__":
    from pysat.formula import CNF

    resolution = 10000
    duration = 5
    file_name = "coloring_basic"
    f1 = CNF(from_file=f"test_files/{file_name}.cnf")
    result_name = f"min_sqr_{file_name}_solveivp_{duration}s_{resolution}res"

    import time

    ode_sat = OdeSat(clauses=f1.clauses, resolution=resolution, time=duration)
    start = time.time()
    sTs, aTs = ode_sat.integrate()
    end = time.time()
    print(sTs)
    last = sTs[-1]
    print(last.tolist())
    print([index + 1 if value > 0 else -index - 1 for index, value in enumerate(last)])
    plt.plot(sTs)
    plt.ylim(-1.1, 1.1)
    plt.legend(np.arange(ode_sat.I) + 1)
    plt.savefig(f"out/{result_name}.png")

    # 17 4 8 14 2 6
    plt.clf()
    plt.plot(sTs[:, [16, 3, 7, 13, 1, 5]])
    plt.legend([17, 4, 8, 14, 2, 6])
    plt.ylim(-1.1, 1.1)
    plt.savefig(f"out/{result_name}_sel.png")

    # 17 4 8 14 2 6
    plt.clf()
    KTs = np.array(KTs)
    plt.plot(KTs[:, [7, 18]])
    plt.legend([8, 19])
    plt.savefig(f"out/{result_name}_K_sel.png")

    plt.clf()
    vs = [ode_sat.V(s, a) for s, a in zip(sTs, aTs)]
    plt.plot(vs)
    plt.savefig(f"out/{result_name}_v.png")

    plt.clf()
    es = [ode_sat.E(s) for s in sTs]
    plt.plot(es)
    plt.savefig(f"out/{result_name}_e.png")

    plt.clf()
    plt.plot(KTs)
    plt.savefig(f"out/{result_name}_K.png")

    plt.clf()
    plt.plot(aTs)
    plt.ylim(0, 10)
    # plt.yscale("log")
    plt.savefig(f"out/{result_name}_a.png")

    print(f"Time: {end - start} s, final V: {vs[-1]}")

    print(verify(f"test_files/{file_name}.cnf", last.tolist()))

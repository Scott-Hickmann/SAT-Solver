import numpy as np
from pysat.formula import CNF
from matplotlib import pyplot as plt
from solver import verify

test_name = "complex"
file_name = f"test_files/{test_name}.cnf"
cnf = CNF(from_file=file_name)

np.random.seed(0)


class Variable:
    BUCKET_SIZE = 1000

    def __init__(self, index):
        self.index = index
        self.bucket_ones = Variable.BUCKET_SIZE // 2

    def select_bit(self):
        bit = 1 if np.random.randint(0, Variable.BUCKET_SIZE) < self.bucket_ones else 0
        self.value = bit

    def add(self, change):
        self.bucket_ones += change
        self.bucket_ones = max(min(self.bucket_ones, Variable.BUCKET_SIZE), 0)

    def __repr__(self):
        return f"V{self.index + 1}={self.bucket_ones / Variable.BUCKET_SIZE}"


class Clause:
    def __init__(self, clause):
        self.clause = clause
        self.a = 0  # represents a bit shift

    def variable_details(self):
        for v in self.clause:
            vindex = abs(v) - 1
            vpositive = v > 0
            yield vindex, vpositive

    def variable_satisfied(self, vindex, vpositive):
        vvalue = variables[vindex].value
        return (vvalue and vpositive) or (not vvalue and not vpositive)

    # True if clause is satisfied
    def satisfied(self, ignore_vindex=None):
        for vindex, vpositive in self.variable_details():
            if ignore_vindex == vindex:
                continue
            if self.variable_satisfied(vindex, vpositive):
                return True
        return False

    def update(self):
        satisfied = self.satisfied()
        K = int(not satisfied)  # 1 if clause unsatisfied with all variables

        dsdt_contributions = np.zeros(cnf.nv, dtype=int)

        for vindex, vpositive in self.variable_details():
            satisfied = self.satisfied(ignore_vindex=vindex)
            Ki = int(not satisfied)  # 1 if clause not satisfied without variable i
            ci = int(1 if vpositive else -1)
            contribution = (K * Ki * ci) << max(round(self.a), 0)
            dsdt_contributions[vindex] = contribution

        # update a
        if K:
            self.a += 1.0/50

        return dsdt_contributions


variables = np.array([Variable(i) for i in range(cnf.nv)])
clauses = np.array([Clause(clause) for clause in cnf.clauses])


def update():
    dsdt_contributions = np.zeros(cnf.nv, dtype=int)
    for run in range(50):
        for variable in variables:
            variable.select_bit()

        for clause in clauses:
            dsdt_contributions += clause.update()

    dsdt_contributions = dsdt_contributions // 50
    for i, contribution in enumerate(dsdt_contributions):
        variables[i].add(contribution)

    max_a = max([clause.a for clause in clauses])
    diff_a = max(max_a - 5, 0)
    for clause in clauses:
        clause.a -= diff_a


history = []
history_a = []
for i in range(1000):
    update()
    history_a.append([c.a for c in clauses])
    history.append([v.bucket_ones / Variable.BUCKET_SIZE for v in variables])

print(variables)

result = [1 if x > 0 else -1 for x in history[-1]]

result_fname = f"out/digital/{test_name}.png"
result_fname_a = f"out/digital/{test_name}_a.png"

plt.plot(history)
plt.legend([variable.index for variable in variables])
plt.savefig(result_fname)

plt.plot(history_a)
plt.ylim(0, 5)
plt.legend([i + 1 for i in range(len(clauses))])
plt.savefig(result_fname_a)

print(f"saved output graph to {result_fname}")

solved, failed_clauses = verify(file_name, result)

print("Valid solution:", solved)
print("Failed clauses:", failed_clauses)

failed_clauses_redundant = []
for index, clause in enumerate(clauses):
    if not clause.satisfied():
        failed_clauses_redundant.append(index + 1)
print("Failed clauses redundant:", failed_clauses_redundant)

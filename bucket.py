import numpy as np
from pysat.formula import CNF
from matplotlib import pyplot as plt
from solver import count_solved


class Bucket:
    def __init__(self, initial_size=500, initial_value=0.5):
        self.size = initial_size
        self.kept_rate = 0.999
        self.bucket = list(
            np.random.choice(
                [0, 1], size=initial_size, p=[1 - initial_value, initial_value]
            )
        )

    def add_bit(self, bit: bool):
        self.bucket.append(bit)

    def get_bit(self) -> bool:
        b = 0
        while len(self.bucket) > self.size:
            b = self.bucket.pop(np.random.randint(len(self.bucket)))
        self.leak()
        return b

    def leak(self):
        self.size *= self.kept_rate
        self.bucket = self.bucket[: int(self.size)]

    def get_value(self) -> float:
        return sum(self.bucket) / len(self.bucket)

    def __repr__(self) -> str:
        return f"Bucket {np.round(self.get_value(), 2)} ({sum(self.bucket)} / {len(self.bucket)})"


def term_unsatisfied(variable_normal: bool, value: bool):
    """
    variableNormal = True if variable is normal else False
    variable = True if variable is True else False
    returns True if that part of the clause is not satisfied
    """
    return variable_normal ^ value


# Returns 1 if the clause is not satisfied
def clause_unsatisfied(terms_unsatisfied: list):
    return np.all(terms_unsatisfied)


# randomly pulls from the varible XOR a clause
def get_new_values(
    clauses_unsatisfied: list, clauses_variables: list, value: bool, variable: int
):
    updated_variable_pool = []
    a_clause_unsatisfied = False
    for clause_variables, clause_unsatisfied in zip(
        clauses_variables, clauses_unsatisfied
    ):
        if clause_unsatisfied:
            a_clause_unsatisfied = True
        if variable in clause_variables:
            updated_variable_pool.append(clause_unsatisfied ^ value)
    return updated_variable_pool, not a_clause_unsatisfied


problem_file = "test_files/coloring_basic.cnf"
problem = CNF(from_file=problem_file)

num_variables = problem.nv
clauses = [[variable > 0 for variable in clause] for clause in problem.clauses]
clauses_variables = [
    [abs(variable) - 1 for variable in clause] for clause in problem.clauses
]

buckets = [Bucket(initial_size=500) for _ in range(num_variables)]


def update():
    values = [bucket.get_bit() for bucket in buckets]
    clauses_unsatisfied = []
    for clause_variables, clause in zip(clauses_variables, clauses):
        terms_unsatisfied = [
            term_unsatisfied(variable_normal, values[variable])
            for variable, variable_normal in zip(clause_variables, clause)
        ]
        clauses_unsatisfied.append(clause_unsatisfied(terms_unsatisfied))
    for variable, value in enumerate(values):
        new_values, done = get_new_values(
            clauses_unsatisfied, clauses_variables, value, variable
        )
        for new_value in new_values:
            buckets[variable].add_bit(new_value)

    return np.any([bucket.size == 1 for bucket in buckets])


history = []
for i in range(100000):
    done = update()
    if done:
        print(f"\nfinished in {i} iterations")
        break
    history.append([bucket.get_value() for bucket in buckets])
    if i % 2000 == 0:
        print([bucket for bucket in buckets])
print([bucket for bucket in buckets])

result = [
    variable + 1 if round(bucket.get_value()) > 0 else -(variable + 1)
    for variable, bucket in enumerate(buckets)
]
print(result)
print(f"clauses correct: {count_solved(problem_file, result)} / {len(problem.clauses)}")

plt.plot(history)
plt.ylim([-0.1, 1.1])
plt.legend(range(1, len(buckets) + 1))
savepath = f"out/{problem_file.split('/')[-1][:-4]}.png"
plt.savefig(savepath)
print("Saved to " + savepath)

"""
if all termSatisfied is False, then we need to "kick" the system
by taking all the variables that fed into that clause, 
and adding the opposite values into those respective variable buckets
"""

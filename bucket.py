import numpy as np
from pysat.formula import CNF


class Bucket:
    def __init__(self, initial_size=1000, initial_value=0.5):
        self.bucket = list(
            np.random.choice(
                [0, 1], size=initial_size, p=[1 - initial_value, initial_value]
            )
        )

    def add_bit(self, bit: bool):
        self.bucket.append(bit)

    def get_bit(self) -> bool:
        return self.bucket.pop(np.random.randint(len(self.bucket)))

    def get_value(self) -> float:
        return sum(self.bucket) / len(self.bucket)

    def __repr__(self) -> str:
        return f"Bucket {self.get_value()} ({sum(self.bucket)} / {len(self.bucket)})"


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
def get_new_value(
    clauses_unsatisfied: list, clauses_variables: list, value: bool, variable: int
):
    updated_variable_pool = []
    for clause_variables, clause_unsatisfied in zip(
        clauses_variables, clauses_unsatisfied
    ):
        if variable in clause_variables:
            updated_variable_pool.append(clause_unsatisfied ^ value)
    return np.random.choice(updated_variable_pool)


problem = CNF(from_file="test_files/complex.cnf")

num_variables = problem.nv
clauses = [[variable > 0 for variable in clause] for clause in problem.clauses]
clauses_variables = [
    [abs(variable) - 1 for variable in clause] for clause in problem.clauses
]

buckets = [Bucket() for _ in range(num_variables)]


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
        new_value = get_new_value(
            clauses_unsatisfied, clauses_variables, value, variable
        )
        buckets[variable].add_bit(new_value)


print([bucket for bucket in buckets])
for i in range(10000):
    update()
    print([bucket for bucket in buckets])

result = [round(bucket.get_value()) for bucket in buckets]
print(result)

"""
if all termSatisfied is False, then we need to "kick" the system
by taking all the variables that fed into that clause, 
and adding the opposite values into those respective variable buckets
"""

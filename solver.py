from pysat.formula import CNF, WCNF
from pysat.examples.fm import FM


def solve(file_name: str):
    f1 = CNF(from_file=file_name)
    wf1 = WCNF()
    wf1.extend(f1.clauses, weights=[1 for _ in range(len(f1.clauses))])
    with FM(wf1) as fm:
        print(fm.compute())
        return fm.model


def count_solved(file_name: str, result: list):
    f1 = CNF(from_file=file_name)
    correct_clauses = 0
    for clause in f1.clauses:
        clause_result = False
        for variable in clause:
            value = result[abs(variable) - 1] > 0
            if variable > 0:
                clause_result = clause_result or value
            else:
                clause_result = clause_result or not value
        correct_clauses += 1 if clause_result else 0
    return correct_clauses


def verify(file_name: str, result: list):
    f1 = CNF(from_file=file_name)
    for clause in f1.clauses:
        clause_result = False
        for variable in clause:
            value = result[abs(variable) - 1] > 0
            if variable > 0:
                clause_result = clause_result or value
            else:
                clause_result = clause_result or not value
        if not clause_result:
            return False
    return True


if __name__ == "__main__":
    print(
        verify(
            "test_files/coloring_basic.cnf",
            [
                -1,
                -2,
                3,
                4,
                -5,
                -6,
                -7,
                -8,
                -9,
                10,
                11,
                -12,
                -13,
                14,
                -15,
                -16,
                17,
                18,
                19,
                -20,
                -21,
                22,
            ],
        )
    )
    # soln = solve("test_files/complex.cnf")

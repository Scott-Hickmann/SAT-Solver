from pysat.formula import CNF, WCNF
from pysat.examples.fm import FM


def solve(file_name: str):
    f1 = CNF(from_file=file_name)
    wf1 = WCNF()
    wf1.extend(f1.clauses, weights=[1 for _ in range(len(f1.clauses))])
    with FM(wf1) as fm:
        print(fm.compute())
        return fm.model


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
    soln = solve("test_files/factor21.cnf")
    print(soln)

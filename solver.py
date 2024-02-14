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
    """
    result should be a list of the variable values,
    with -1 for false and 1 for true.
    e.g. x1 = True, x2 = False, x3 = True, x4 = False
    result = [1, -1, 1, -1]
    """
    f1 = CNF(from_file=file_name)
    failed_clauses = []
    for index, clause in enumerate(f1.clauses):
        clause_result = False
        for variable in clause:
            value = result[abs(variable) - 1] > 0
            if variable > 0:
                clause_result = clause_result or value
            else:
                clause_result = clause_result or not value
        if not clause_result:
            failed_clauses.append(index + 1)
    if len(failed_clauses) > 0:
        return False, failed_clauses
    return True, failed_clauses


if __name__ == "__main__":
    soln = solve("test_files/complex.cnf")
    print(soln)

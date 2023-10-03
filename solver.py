from pysat.formula import CNF, WCNF
from pysat.examples.fm import FM

def solve_max(file_name: str):
    f1 = CNF(from_file=file_name)
    wf1 = WCNF()
    wf1.extend(f1.clauses, weights=[1 for _ in range(len(f1.clauses))])
    with FM(wf1) as fm:
        print(fm.compute())
        return fm.model

if __name__ == "__main__":
    soln = solve_max('test_files/simple.cnf')
    print(soln)
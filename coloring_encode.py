from pysat.formula import CNF
import numpy as np

from coloring_utils import graph, get_color_variables, get_pair_variable

# Variable a_i_j, b_i_j

# Use np.where to find the row and column indices where the matrix has a value of 1
rows, cols = np.where(graph == 1)

# Create the pairs of neighbors
pairs = np.array([[r, c] for r, c in zip(rows, cols) if r < c])

print("Neighbors", [[chr(65 + number) for number in pair] for pair in pairs])

clauses = np.empty((4 * len(pairs), 3), dtype=int)
for pair_index, pair in enumerate(pairs):
    a_i, b_i = get_color_variables(pair[0])
    a_j, b_j = get_color_variables(pair[1])
    c = get_pair_variable(pair_index)
    clause_index = 4 * pair_index
    print(clause_index)
    clauses[clause_index] = [c, a_i, a_j]
    clauses[clause_index + 1] = [c, -a_i, -a_j]
    clauses[clause_index + 2] = [-c, b_i, b_j]
    clauses[clause_index + 3] = [-c, -b_i, -b_j]
    # print("Pair", pair_index, "Variables", a_i, b_i, a_j, b_j, c)

print("Clauses", clauses)

f1 = CNF()
f1.extend(clauses)
f1.to_file('test_files/coloring.cnf')
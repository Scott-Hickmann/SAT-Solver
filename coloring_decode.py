import numpy as np
from pysat.formula import CNF

from coloring_utils import graph, parse_solution_variable

solution = [-1, 2, -3, -4, 5, -6, -7, 8, -9, -10, 11, -12, 13, -14, 15, -16, 17, -18, -19, -20, -21, -22]

color_sets = [set() for _ in range(4)]

for node in range(graph.shape[0]):
    node, color_index = parse_solution_variable(node, solution)
    color_sets[color_index].add(node)

print("Color sets", color_sets)
import numpy as np

# Adjascent matrix
graph = np.array([
    [0, 1, 1, 0, 1, 0],
    [1, 0, 1, 1, 0, 1],
    [1, 1, 0, 1, 1, 0],
    [0, 1, 1, 0, 0, 1],
    [1, 0, 1, 0, 0, 1],
    [0, 1, 0, 1, 1, 0]
])

def get_color_variables(number: int):
    return 2 * number + 1, 2 * number + 2

def get_pair_variable(number: int):
    return 2 * len(graph) + number + 1

def parse_solution_variable(node: int, solution: list[int]):
    a_i_value = 1 if solution[2 * node] > 0 else 0
    b_i_value = 1 if solution[2 * node + 1] > 0 else 0
    color_index = a_i_value + 2 * b_i_value
    return node, color_index
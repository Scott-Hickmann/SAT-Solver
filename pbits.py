import numpy as np
from pysat.formula import CNF

probA = 0.5
probB = 0.5
probC = 0.5
signalLength = 10

file_name = 'test_files/simple.cnf'
cnf = CNF(from_file=file_name)
clauses = cnf.clauses

signals = np.array([
    np.array([1 if np.random.random() < probA else 0 for i in range(signalLength)]),
    np.array([1 if np.random.random() < probB else 0 for i in range(signalLength)]),
    np.array([1 if np.random.random() < probC else 0 for i in range(signalLength)])
])
print(signals[0])
print(signals[1])
print(signals[2])

clause_signals = []
for bits in zip(*signals):
    clause_signal = []
    for clause in clauses:
        maybe_negated_values = []
        for variable in clause:
            variable_index = abs(variable) - 1
            if variable > 0:
                maybe_negated = 1 - bits[variable_index]
            else:
                maybe_negated = bits[variable_index]
            maybe_negated_values.append(maybe_negated)
        clause_output = np.min(maybe_negated_values) # AND
        clause_signal.append(clause_output)
    clause_signals.append(clause_signal)
    # output_signal.append(np.max(clause_outputs)) # OR

print(clause_signals)

    # print(Ki)
    # Ks.append(Ki)
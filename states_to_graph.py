import numpy as np

# List of state pairs
state_pairs = [
    "AK,WA", "AL,FL", "AL,GA", "AL,MS", "AL,TN", "AR,LA", "AR,MO", "AR,MS", 
    "AR,OK", "AR,TN", "AR,TX", "AZ,CA", "AZ,CO", "AZ,NM", "AZ,NV", "AZ,UT", 
    "CA,HI", "CA,NV", "CA,OR", "CO,KS", "CO,NE", "CO,NM", "CO,OK", "CO,UT", 
    "CO,WY", "CT,MA", "CT,NY", "CT,RI", "DC,MD", "DC,VA", "DE,MD", "DE,NJ", 
    "DE,PA", "FL,GA", "GA,NC", "GA,SC", "GA,TN", "IA,IL", "IA,MN", "IA,MO", 
    "IA,NE", "IA,SD", "IA,WI", "ID,MT", "ID,NV", "ID,OR", "ID,UT", "ID,WA", 
    "ID,WY", "IL,IN", "IL,KY", "IL,MO", "IL,WI", "IN,KY", "IN,MI", "IN,OH", 
    "KS,MO", "KS,NE", "KS,OK", "KY,MO", "KY,OH", "KY,TN", "KY,VA", "KY,WV", 
    "LA,MS", "LA,TX", "MA,NH", "MA,NY", "MA,RI", "MA,VT", "MD,PA", "MD,VA", 
    "MD,WV", "ME,NH", "MI,OH", "MI,WI", "MN,ND", "MN,SD", "MN,WI", "MO,NE", 
    "MO,OK", "MO,TN", "MS,TN", "MT,ND", "MT,SD", "MT,WY", "NC,SC", "NC,TN", 
    "NC,VA", "ND,SD", "NE,SD", "NE,WY", "NH,VT", "NJ,NY", "NJ,PA", "NM,OK", 
    "NM,TX", "NM,UT", "NV,OR", "NV,UT", "NY,PA", "NY,VT", "OH,PA", "OH,WV", 
    "OK,TX", "OR,WA", "PA,WV", "SD,WY", "TN,VA", "UT,WY", "VA,WV"
]

# Extract unique states and sort them
states = sorted(list(set(state for pair in state_pairs for state in pair.split(','))))

# Create a mapping from state name to index
state_to_index = {state: index for index, state in enumerate(states)}

# Initialize adjacency matrix with zeros
n = len(states)
adj_matrix = np.zeros((n, n), dtype=int)

# Fill in the adjacency matrix
for pair in state_pairs:
    state1, state2 = pair.split(',')
    i, j = state_to_index[state1], state_to_index[state2]
    adj_matrix[i, j] = 1
    adj_matrix[j, i] = 1  # because it's an undirected graph

print(state_to_index)
# print(adj_matrix.tolist())
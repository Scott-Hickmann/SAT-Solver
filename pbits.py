import numpy as np
from matplotlib import pyplot as plt

# returns a vector of specified size
# where each element has p probability of being a 1
# def ber_stream(size, p):
#     return np.random.binomial(size=size, n=1, p=p)


# p probability of returning a True, else False
def ber(p):
    return p > np.random.uniform(0, 1)


def minus(a, b):
    return a or not b


def add(a, b):
    return np.random.choice(1, (a, b))


def mul(a, b):
    return a or b


def moving_average(x, w):
    return np.convolve(x, np.ones(w), "valid") / w


# TODO implement ode dy/dx = target - y
STREAM_LENGTH = 100000
target = 0.8
x = 0
y = 0
history = []
history_x = []
window_size = 1000  # Window size for rolling average

# x(0) = 1
# y(0) = 1
# x' = - x - y + 0.6
# y' = - x - y + 0.6
# --> 0.3, 0.3

0.6111

0000
0101

0001

# "1 + 0"  -> 0 + 0.5 -> 0.25 -> 0.5 but should be 1

1110 (-0.5)


x_0 = 00001000000000000000000000000000000

x_1 = 0010101010100101

x_2 = 110101010010101

and (x_0 \/ x_1 \/ x_2) = 0000000000010000001000000000

for _ in range(STREAM_LENGTH):
    dx = minus(minus(ber(0.6), x), y)
    dy = minus(minus(ber(0.6), x), y)
    x = add(x, dx) if ber(0.01) else x
    y = add(y, dy) if ber(0.01) else y
    history_x.append(x)
    history.append(y)

# Convert boolean values to integers and calculate rolling average
history_int = np.array(history, dtype=int)
rolling_avg = np.convolve(history_int, np.ones(window_size) / window_size, mode="valid")
rolling_avg_x = np.convolve(
    np.array(history_x, dtype=int), np.ones(window_size) / window_size, mode="valid"
)

print(history)

# Plotting
plt.plot(rolling_avg)
plt.plot(rolling_avg_x)
plt.title("Rolling Average of History")
plt.xlabel("Iteration")
plt.ylabel("Average Value")
plt.ylim(-0.1, 1.1)
plt.savefig("fig.jpg")

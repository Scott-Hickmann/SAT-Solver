from sympy import (
    symbols,
    Sum,
    Product,
    diff,
    IndexedBase,
    Idx,
    piecewise_fold,
    latex,
    Min,
)

# Redefine symbols with indexed bases
M, N = symbols("M N", integer=True, positive=True)
m, j = symbols("m j", integer=True, positive=True)
a = IndexedBase("a", shape=(M,))
c = IndexedBase("c", shape=(M, N))
s = IndexedBase("s", shape=(N,))
i = Idx("i")
m = Idx("m")
j = Idx("j")

# K^2
# K_m = Product(1 - c[m, j] * s[j], (j, 1, N))
# V = Sum(a[m] * K_m**2, (m, 1, M))

# K
K_m = Product(1 - c[m, j] * s[j], (j, 1, N))
V = Sum(a[m] * K_m, (m, 1, M))

# K has min
# K_m = Min(1 - c[m, j] * s[j], (j, 1, N))
# K_m = Min(*[1 - c[m, j] * s[j] for j in range(1, N + 1)])
# V = Sum(a[m] * K_m, (m, 1, M))

# Compute derivative of V with respect to s_i
dV_dsi = -diff(V, s[i])

result = piecewise_fold(dV_dsi.simplify())

print(latex(result))

import pickle
import sympy as sp

A1,A2,b12,J12 = sp.symbols('A_1 A_2 b_{12} J_{12}',real=True)
t = sp.symbols('t',real=True,positive=True)

expr = 4*((b12)**2*(A1-A2)**2)/((b12)**2 + (A1-A2)**2)**2*(sp.sin(t*sp.Rational(1,4)*sp.sqrt(A1**2-2*A1*A2+A2**2+(b12)**2)))**4
expr_j = 4*((b12-2*J12)**2*(A1-A2)**2)/((b12-2*J12)**2 + (A1-A2)**2)**2*(sp.sin(t*sp.Rational(1,4)*sp.sqrt(A1**2-2*A1*A2+A2**2+(b12-2*J12)**2)))**4

with open("I1o2I1o2WtTCL2.txt", "wb") as f:
    pickle.dump(expr, f, protocol=pickle.HIGHEST_PROTOCOL)

with open("I1o2I1o2WtTCL2J.txt", "wb") as f:
    pickle.dump(expr_j, f, protocol=pickle.HIGHEST_PROTOCOL)
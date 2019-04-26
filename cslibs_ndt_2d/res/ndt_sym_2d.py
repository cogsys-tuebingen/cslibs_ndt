#!/usr/bin/python3
from sympy import *

tx,ty,gamma = symbols('tx ty gamma')
x,y = symbols('x y')

R_z = Matrix([[cos(gamma),-sin(gamma)],[sin(gamma), cos(gamma)]])

v = Matrix([[x],[y]])
t = Matrix([[tx],[ty]])

v_prime = (R_z * v) + t

d = Matrix([tx,ty,gamma])
J = v_prime.jacobian(d)

init_printing()

pprint(v_prime)
pprint(J)

print(latex(v_prime[0]))
print("-------------------------")
print(latex(v_prime[1]))
print("-------------------------")
print(latex(J.col(0)))
#print("-------------------------")
#print(latex(J.jacobian(d)))

hessian = []
for i in range(3):
    print("-------------------------")
    print(latex(J.col(i).jacobian(d)))
    hessian.append(J.col(i).jacobian(d))
pprint(hessian)

#!/usr/bin/python3
from sympy import *
from sympy.utilities.codegen import codegen


tx,ty,tz,alpha,beta,gamma = symbols('tx ty tz alpha beta gamma')
x,y,z = symbols('x y z')

R_x = Matrix([[1, 0, 0],[0, cos(alpha), -sin(alpha)],[0,sin(alpha),cos(alpha)]])
R_y = Matrix([[cos(beta), 0, sin(beta)],[0, 1, 0],[-sin(beta),0,cos(beta)]])
R_z = Matrix([[cos(gamma),-sin(gamma), 0],[sin(gamma), cos(gamma), 0],[0,0,1]])

v = Matrix([[x],[y],[z]])
t = Matrix([[tx],[ty],[tz]])

R = R_z * R_y * R_x

v_prime = (R * v) + t

d = Matrix([tx,ty,tz,alpha,beta,gamma])
J = v_prime.jacobian(d)

init_printing()

pprint(v_prime)
pprint(J)

print(latex(v_prime[0]))
print("-------------------------")
print(latex(v_prime[1]))
print("-------------------------")
print(latex(v_prime[2]))
print("-------------------------")
print(latex(J.col(0)))
print("-------------------------")
print(latex(J.col(1)))
print("-------------------------")
print(latex(J.col(2)))
hessian = []
for i in range(6):
    print("-------------------------")
    print(latex(J.col(i).jacobian(d)))
    hessian.append(J.col(i).jacobian(d))

print("output is row major")
[(c_name, c_code), (h_name, c_header)] = codegen([('jacobian', J.col(4))], 'C')
print(c_code)

[(c_name, c_code), (h_name, c_header)] = codegen([('hessian', hessian[5].col(5))], 'C')
print(c_code)

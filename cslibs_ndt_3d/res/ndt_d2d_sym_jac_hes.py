#!/usr/bin/python3
from sympy import *
from sympy.utilities.codegen import codegen

def mat_diff_by_scalar(M, a):
	s = M.shape
	D = zeros(s[0],s[1])

	for i in range(s[0]):
		for j in range(s[1]):
			D[i,j] = diff(M[i,j], a)

	return D



tx,ty,tz,alpha,beta,gamma = symbols('tx ty tz alpha beta gamma')
x,y,z = symbols('x y z')

# ROTATION MATRICES / the general rotation matrix
R_x = Matrix([[1, 0, 0],[0, cos(alpha), -sin(alpha)],[0,sin(alpha),cos(alpha)]])
R_y = Matrix([[cos(beta), 0, sin(beta)],[0, 1, 0],[-sin(beta),0,cos(beta)]])
R_z = Matrix([[cos(gamma),-sin(gamma), 0],[sin(gamma), cos(gamma), 0],[0,0,1]])
R = R_z * R_y * R_x
C = Matrix([[x],[y],[z]]) * Matrix([[x,y,z]])

# THE TRANSLATION
t = Matrix([[tx],[ty],[tz]])

# THE INPUT VECTOR AND ITS TRANSFORMED VERSION
v = Matrix([[x],[y],[z]])
v_prime = (R * v) + t
C_prime = (R.transpose() * C * R)

# TO DERIVE AFTER THESE COMPONENTS
d = Matrix([tx,ty,tz,alpha,beta,gamma])
# DERIVING ONLY v_prime IS SUFFICIENT, SINCE THE MEAN OF THE COUNTER DISTRIBUTION IS NOT
# MOVED AND WILL BECOME 0 => AS A RESULT, THE JACOBIAN IS THE SAME AS FOR THE SAMPLE BASED VERSION
J = v_prime.jacobian(d)
J_mat = []

init_printing()



for var in [tx,ty,tz,alpha,beta,gamma]:
	Z = mat_diff_by_scalar(C_prime, var)
	R_prime = mat_diff_by_scalar(R, var)
	Z_prime = R_prime.transpose() * C * R + R.transpose() * C * R_prime

	if Z_prime == Z:
		print("OK")

	J_mat.append(R_prime)
	pprint(R_prime)

for i in range(3,6):
	[(c_name, c_code), (h_name, c_header)] = codegen([('jacobian', J_mat[i])], 'C')
	print(c_code)



exit(0)

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

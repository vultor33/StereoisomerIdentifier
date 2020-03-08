import sympy
import pandas as pd

a, b, c, d, e, f, g, h, i, j = sympy.symbols('a b c d e f g h i j')
sympy.init_printing(use_unicode=True)
def expansao_polinomios(indice_ciclos):
    return sympy.Poly(indice_ciclos).as_expr()
    
# Definicao das variaves que representam as 10 cores possiveis
x1 = a + b + c + d + e + f + g + h + i + j
x2 = a**2 + b**2 + c**2 + d**2 + e**2 + f**2 + g**2 + h**2 + i**2 + j**2
x3 = a**3 + b**3 + c**3 + d**3 + e**3 + f**3 + g**3 + h**3 + i**3 + j**3
x4 = a**4 + b**4 + c**4 + d**4 + e**4 + f**4 + g**4 + h**4 + i**4 + j**4
x5 = a**5 + b**5 + c**5 + d**5 + e**5 + f**5 + g**5 + h**5 + i**5 + j**5
x6 = a**6 + b**6 + c**6 + d**6 + e**6 + f**6 + g**6 + h**6 + i**6 + j**6
x7 = a**7 + b**7 + c**7 + d**7 + e**7 + f**7 + g**7 + h**7 + i**7 + j**7
x8 = a**8 + b**8 + c**8 + d**8 + e**8 + f**8 + g**8 + h**8 + i**8 + j**8
x9 = a**9 + b**9 + c**9 + d**9 + e**9 + f**9 + g**9 + h**9 + i**9 + j**9
x10 = a**10 + b**10 + c**10 + d**10 + e**10 + f**10 + g**10 + h**10 + i**10 + j**10

# Definicao dos coeficientes desejados
coeficientes_6 =  [a**6, 
                a**5*b,
                a**4*b**2,
                a**4*b*c,
                a**3*b**3,
                a**3*b**2*c,
                a**3*b*c*d,
                a**2*b**2*c**2,
                a**2*b**2*c*d,
                a**2*b*c*d*e,
                a*b*c*d*e*f]
    
coeficientes_10 = [a**10,
                   a**9*b,
                   a**8*b**2,
                   a**8*b*c,
                   a**7*b**3,
                   a**7*b**2*c,
                   a**7*b*c*d,
                   a**6*b**4,
                   a**6*b**3*c,
                   a**6*b**2*c**2,
                   a**6*b**2*c*d,
                   a**6*b*c*d*e,
                   a**5*b**5,
                   a**5*b**4*c,
                   a**5*b**3*c**2,
                   a**5*b**3*c*d,
                   a**5*b**2*c**2*d,
                   a**5*b**2*c*d*e,
                   a**5*b*c*d*e*f,
                   a**4*b**4*c**2,
                   a**4*b**4*c*d,
                   a**4*b**3*c**3,
                   a**4*b**3*c**2*d,
                   a**4*b**3*c*d*e,
                   a**4*b**2*c**2*d**2,
                   a**4*b**2*c**2*d*e,
                   a**4*b**2*c*d*e*f,
                   a**4*b*c*d*e*f*g,
                   a**3*b**3*c**3*d,
                   a**3*b**3*c**2*d**2,
                   a**3*b**3*c**2*d*e,
                   a**3*b**3*c*d*e*f,
                   a**3*b**2*c**2*d**2*e,
                   a**3*b**2*c**2*d*e*f,
                   a**3*b**2*c*d*e*f*g,
                   a**3*b*c*d*e*f*g*h,
                   a**2*b**2*c**2*d**2*e**2,
                   a**2*b**2*c**2*d**2*e*f,
                   a**2*b**2*c**2*d*e*f*g,
                   a**2*b**2*c*d*e*f*g*h,
                   a**2*b*c*d*e*f*g*h*i,
                   a*b*c*d*e*f*g*h*i*j]




    
# TESTE

# PPY-6
#indice_ciclos_rot = (x1**6 + x1*x5 + x1*x5 + x1*x5 + x1*x5)
#indice_ciclos_completo = (indice_ciclos_rot + x1*x1*x2*x2 
#                          + x1*x2*x1*x2 + x1*x2*x2*x1 + x1*x2*x2*x1 
#                          + x1*x2*x2*x1)
#coeficientes = coeficientes_6
#termo_rot = 1/5
#termo_comp = 1/10

# RESULTADO
#0formula;1total;2chiral;3achiral
#a**6;1;0;1
#a**5*b;2;0;2
#a**4*b**2;3;0;3
#a**4*b*c;6;4;2
#a**3*b**3;4;0;4
#a**3*b**2*c;12;8;4
#a**3*b*c*d;24;24;0
#a**2*b**2*c**2;18;12;6
#a**2*b**2*c*d;36;32;4
#a**2*b*c*d*e;72;72;0
#a*b*c*d*e*f;144;144;0
from sympy import Symbol, Derivative
import numpy as np

#
#
#  Trabalho2:Parte 3 -Encontra uma raiz aproximada para os sistemas de equações abaixo,usando o metodo de Newton para sistemas não lineares.
#
#  F(x1,x2) = 6x1^3 + x1*x2 - 3x2^3 - 4 = 0
#  F(x1,x2) = x2^2 - 18x1*x2^2 + 16x2^3 + 1 = 0 
#
#  
#
#
#
#
#  Autor       : Marlon Brendo Ramos
#  Instituição : Universidade Federal de Uberlândia
#  Curso       : Sistemas de Informação
#  Disciplina  : Métodos Numéricos
#  Tecnologias : Python3.7
#
#  --------------------------------------------------------------------------------------------------------------------------
#  Esse programa permite a resolução do sistema de equações não lineares abaixo,usando o método de Newton
#
#  F(x1,x2) = 6x1^3 + x1*x2 - 3x2^3 - 4 = 0
#  F(x1,x2) = x2^2 - 18x1*x2^2 + 16x2^3 + 1 = 0
#
#  O programa imprime todas as iterações na saída padrão e seus respectivos valores de k , x^k+1 e || x^k+1 - x^k || ,
#  considerando a precisão de e=10^-3 ou 0.001.
#
#
#
#  Exemplo de saída:
#
#--------Iteracao k ------------- 0
#
#   Xk =  [2, 2]
#   X^k+1 =  [1.36864407 1.33474576]
#   || x^k+1 - x^k || =  0.6652542372881356
#
#----------------------------------
#...
#
#

#Calcula  a matriz jacobiana, utilizando 
#derivadas parciais em relaçao a x e y e monta
# uma matriz com as derivadas
def calculaMatrizjacobiana():
    
    f1 = '6*x**3 + x*y - 3*y**3 -4'
    f2 = 'x**2 - 18*x*y**2 + 16*y**3 +1'
    x1= Symbol('x')
    x2= Symbol('y')
    
    matrizjacobiana = np.zeros((2,2))
    matrizjacobiana=[

                      [ str( Derivative.doit(Derivative(f1, x1))),   str( Derivative.doit(Derivative(f1, x2)) ) ],
                      [ str( Derivative.doit(Derivative(f2, x1)) ) , str( Derivative.doit(Derivative(f2, x2)) ) ]

    ]

    return matrizjacobiana

#Define J(x^k) iniciada por valores de  xk
#tal que x0[2,2]
def defineJxk(xk,matrizjacobiana):
    
    x=xk[0]
    y=xk[1]
    Fxk=[["",""],["",""]]
    for i in range(2):
        for j in range(2):
                Fxk[i][j]=eval(matrizjacobiana[i][j]) 

    return Fxk



#Define -F(x^k).Calcula-se o valor  das funções no xk 
def defineFxk(xk):            
    return [
            ( 6*xk[0]**3 + xk[0]*xk[1] - 3*xk[1]**3 -4)      * -1,
            (xk[0]**2 - 18*xk[0]*xk[1]**2 + 16*xk[1]**3 +1 ) * -1
    ] 




if "__main__" == __name__:

    norma=1
    k=0
    xk=[2,2] #inicia os vetores x0 como [2,2]

    matrizjacobiana=calculaMatrizjacobiana()
   
    while norma > 0.001:
           
        Jxk=defineJxk(xk,matrizjacobiana)
        Fxk=defineFxk(xk)

        #Buscando a inversa das matrizes para achar as soluções
        inv=np.linalg.inv(Jxk) #Aplica a inversa da matriz jacobiana
        sk=np.dot(inv,Fxk) #Encontra as soluções x1,x2 para o sistema 
        xkprox=xk + sk #Soma os vetores de X^k + S^k    
        norma=max(abs(xkprox-xk)) #Encontra a norma fazendo || x^k+1 - x^k ||

        print("--------Iteracao k = {:d}-------------".format(k))
        print("\nXk = ",xk)
        print("X^k+1 = ",xkprox)
        print("|| x^k+1 - x^k || = ",norma)
        print("\n")
        xk=xkprox #atribui o vetor anterior como sendo o próximo
        k+=1



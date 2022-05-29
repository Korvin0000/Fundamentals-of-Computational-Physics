import pylab as pl
import numpy as np

a = 2 #width of well
U_0 = 20 #depth of well
C = np.sqrt(1/(2 * U_0 * a**2))
delta = 10 ** (-9) #accuracy
N = 0

def f(x): #Wavefunction
    
    return np.cos(x) - C * x

def derivative_f(x):
    
    return -np.sin(x) - C

def dichotomy(left, right):
    
    while (right - left > delta):
        middle = (right + left) / 2
        if (f(middle) * f(left) <= 0):
            right = middle
        else:
            left = middle   
    return middle

X = dichotomy(0, np.pi/2)
print('X = ', X)
print('E =', (X**2)/2/a**2 - U_0, "- Dichotomy")      


A = []

for i in np.arange(0, np.pi/2, 0.01):
    
    A.append(abs(derivative_f(i)))
    
G = min(A)
Lambda = 1 / max(A)

def simple_iterations(x_0):
    
    while(abs(f(x_0)) > delta * G):
        
        x_n = x_0 + f(x_0) * Lambda #Lambda < 0 
        x_0 = x_n
        
    return x_n

X = simple_iterations(np.pi/2)
print('X = ', X)
print('E =', (X**2)/2/a**2 - U_0, "- Simple iterations" )



def newton_method(x_0):
    
    while(abs(f(x_0)) > delta * G):
        
        x_n = x_0 - f(x_0) / derivative_f(x_0)
        x_0 = x_n
        
    return x_n

X = newton_method(np.pi/2)
print('X = ', X)
print('E =', (X**2)/2/a**2 - U_0, "- Newton")



P=[]   #array of well depths

for i in np.logspace(-6, 6, 500):
    P.append(i)

def F(x, c):
    return np.cos(x) - c * x

def derivative_F(x, c):
    return -np.sin(x) - c

def dichotomy_with_N(left, right, c):
    N = 0
    while (right - left > delta):
        middle = (right + left) / 2
        if (F(middle, c) * F(left, c) <= 0):
            right = middle
        else:
            left = middle
        N += 1
    return N

def simple_iterations_with_N(x_0, c, G, Lambda):
    N = 0
    while(abs(F(x_0, c)) > delta * G):
        x_n = x_0 + F(x_0, c) * Lambda
        x_0 = x_n
        N += 1
    return N

def newton_method_with_N(x_0, c, G):
    N = 0
    while(abs(F(x_0, c)) > delta * G):
        x_n = x_0 - F(x_0, c) / derivative_F(x_0, c)
        x_0 = x_n
        N += 1
    return N

count_of_iter_dich = []
count_of_iter_sit = []
count_of_iter_newt = []


for i in P:
    
    A = [] #Модули производных при с_i от 0 до pi/2
    for j in np.arange(0, np.pi/2, 0.01):    
        A.append(abs(derivative_F(j, 1/(2 * i * a**2))))
        
    G = min(A) #Min of abs(der)
    L = 1 / max(A)
    count_of_iter_dich.append(dichotomy_with_N(0, np.pi/2, 1/(2 * i * a**2)))
    count_of_iter_sit.append(simple_iterations_with_N(np.pi/2, 1/(2 * i * a**2), G, L))
    count_of_iter_newt.append(newton_method_with_N(np.pi/2, 1/(2 * i * a**2), G))
    

  
pl.figure()
pl.subplot(111)
pl.grid()
pl.semilogx(P, count_of_iter_dich, label = 'Dichotomy', color = 'red')
pl.semilogx(P, count_of_iter_sit, label = 'Simple iterations', color = 'green')
pl.semilogx(P, count_of_iter_newt, label = 'Newton method', color = 'blue')
pl.legend(loc = 'upper right', fontsize = 20)
pl.title('Number of iterations', fontsize = 20)
pl.xlabel('Ln(U)', fontsize = 20)
pl.ylabel('N', fontsize = 20)
pl.show()

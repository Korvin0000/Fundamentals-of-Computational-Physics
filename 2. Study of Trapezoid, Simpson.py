
import numpy as np
import pylab as pl

def f1(x):
    return 1/(1+x*x)

def f2(x):
    return np.power(x, 1/3)*np.exp(np.sin(x))

# Trapezoids

def lab3_1(left, right, N, f):
        s = 0
        dx = (right-left)/N
        x = left
        y1 = f(x)
        for k in range(0, N):
            y2 = f(x+dx)
            s += (y1+y2)*dx*0.5
            y1 = y2
            x += dx
        return s
    
print('X_trapez_1 = ', lab3_1(-1, 1, 512, f1))
print('X_trapez_2 = ', lab3_1(0,1,512,f2))

# Simpson

def lab3_2(left, right, N, f):
    s = 0
    dx = (right - left) / N
    for x in np.arange(left + dx / 2, right - dx / 2 + dx, dx):
        s += dx / 6 * (f(x - dx / 2) + 4 * f(x) + f(x + dx / 2))
    return s

print('X_simpson_1 = ', lab3_2(-1, 1, 512, f1))
print('X_simpson_2 = ', lab3_2(0,1,512,f2))


A = []
for i in range(2, 10+1):
    A.append(2**i)

p = [] #Array of errors
for i in A:
    epsilon_h = abs(np.pi/2 - lab3_1(-1, 1, i, f1))
    epsilon_h_2 = abs(np.pi/2 - lab3_1(-1, 1, 2*i, f1))
    p.append(np.log2(epsilon_h/epsilon_h_2))
    
pl.figure("Trapezoids: 1st")
pl.grid()
pl.plot(A, p, "o-", color = 'red')
pl.title('Dependency of the error on the number of intervals')
pl.xlabel('N')
pl.ylabel('P')
pl.show()
# print(A)


p = []
for i in A:
    epsilon_up = abs(lab3_1(0.01, 1, i ,f2) - lab3_1(0.01, 1, 2*i, f2))
    epsilon_down = abs(lab3_1(0.01, 1, 2*i, f2) - lab3_1(0.01, 1, 4*i, f2))
    p.append(np.log2(epsilon_up/epsilon_down))
    
pl.figure("Trapezoids: 2nd")
pl.grid()
pl.plot(A, p, "o-",  color = 'red')
pl.title('Ошибка от количества интервалов')
pl.xlabel('N')
pl.ylabel('P')
pl.show()


epsilon_h = []
epsilon_h_2 = []

p = []
for i in A:
    epsilo_h = abs(np.pi/2 - lab3_2(-1, 1, i, f1))
    epsilo_h_2 = abs(np.pi/2 - lab3_2(-1, 1, 2*i, f1))
    p.append(np.log2(epsilo_h/epsilo_h_2))
    epsilon_h.append(abs(np.pi/2 - lab3_2(-1, 1, i, f1)))
    epsilon_h_2.append(abs(np.pi/2 - lab3_2(-1, 1, 2*i, f1)))
    
pl.figure("Simpson: 1st")  
pl.grid()
pl.plot(A, p, "o-", color = 'red')
pl.title('Ошибка от количества интервалов')
pl.xlabel('N')
pl.ylabel('P')
pl.show()


p = []
for i in A:
    epsilon_up = abs(lab3_2(0.01, 1, i ,f2) - lab3_2(0.01, 1, 2*i, f2))
    epsilon_down = abs(lab3_2(0.01, 1, 2*i, f2) - lab3_2(0.01, 1, 4*i, f2))
    p.append(np.log2(epsilon_up/epsilon_down))
    
pl.figure("Simpson: 2nd")
pl.grid()
pl.plot(A, p, "o-", color = 'red')
pl.title('Ошибка от количества интервалов')
pl.xlabel('N')
pl.ylabel('P')
pl.show()


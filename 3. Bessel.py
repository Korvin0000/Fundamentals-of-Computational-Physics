import numpy as np
import pylab as pl

def Bessel_0(x, t):
    
    return (1/np.pi)*(np.cos(x*np.sin(t)))

def Bessel_1(x, t):
    
    return (1/np.pi)*np.cos(t-x*np.sin(t))

def Integr_Simpson(x, N, f):
    
    s = 0
    dt = np.pi / N
    for t in np.arange(0, np.pi + dt, dt):
        s += dt / 6 * (f(x, t - dt / 2) + 4 * f(x, t) + f(x, t + dt / 2))
    return s


def dBessel_0 (x, N):
    
        return (Integr_Simpson(x+2*np.pi/N, M, Bessel_0) - Integr_Simpson(x-2*np.pi/N, M, Bessel_0))*N/4/np.pi



Max = []
N = 120
M = 120
Abs_Sum = []
    
sum_0 = abs((-3*Integr_Simpson(0, M, Bessel_0)+4*Integr_Simpson(0+2*np.pi/N, M, Bessel_0) - Integr_Simpson(0+4*np.pi/N, M, Bessel_0))*N/4/np.pi + Integr_Simpson(0, M, Bessel_1))
sum_2_pi = abs((3*Integr_Simpson(2*np.pi, M, Bessel_0) - 4*Integr_Simpson(2*np.pi - 2*np.pi/N, M, Bessel_0) + Integr_Simpson(2*np.pi - 4*np.pi/N, M, Bessel_0))*N/4/np.pi + Integr_Simpson(2*np.pi, M, Bessel_1)) 
    
Abs_Sum.append(sum_0)
Abs_Sum.append(sum_2_pi)
    
for x in np.arange(1, N, 1)*2*np.pi/N:
    sum = abs(dBessel_0(x, N) + Integr_Simpson(x, M, Bessel_1))
    Abs_Sum.append(sum)
#     print(Abs_Sum)
Max.append(max(Abs_Sum))

print(Max)

X = [2**i for i in range(12,13)]

for N in X:

    Abs_Sum = []
    
    sum_0 = abs((-3*Integr_Simpson(0, N, Bessel_0)+4*Integr_Simpson(0+np.pi/N, N, Bessel_0) - Integr_Simpson(0+2*np.pi/N, N, Bessel_0))*N/2/np.pi + Integr_Simpson(0, N, Bessel_1))
    sum_2_pi = abs((3*Integr_Simpson(2*np.pi, N, Bessel_0) - 4*Integr_Simpson(2*np.pi - np.pi/N, N, Bessel_0) + Integr_Simpson(2*np.pi - 2*np.pi/N, N, Bessel_0))*N/2/np.pi + Integr_Simpson(2*np.pi, N, Bessel_1)) 
    
    Abs_Sum.append(sum_0)
    Abs_Sum.append(sum_2_pi)
    
    for x in np.arange(1, 2*N, 1)*np.pi/N:
        sum = abs(dBessel_0(x, N) + Integr_Simpson(x, N, Bessel_1))
      
        Abs_Sum.append(sum)
#     print(Abs_Sum)
    Max.append(max(Abs_Sum))
# print(Max)
# print(X)
pl.figure()
pl.subplot(111)
pl.grid()
pl.semilogy(X, Max, color = 'red')
pl.title('Максимум модуля суммы от количества разбиений')
pl.xlabel('N')
pl.ylabel('Max |Sum|')
pl.show()




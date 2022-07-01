import numpy as np
from math import exp, cos, sin, pow
import matplotlib.pyplot as plt
from scipy.integrate import odeint 
def test(y, t):
    return np.array([1/y[1], -1/y[0]])
def polet(y, t, m, C, rho, S, g, V0, tetha0):
    return np.array([y[2]*cos(y[3]), y[2]*sin(y[3]), (-C*rho*S*y[2]**2)/(2*m) - g*sin(y[3]), -g*cos(y[3])/y[2]])
m=45
C=0.25
rho=1.29
S=0.35
g=9.81
V0=60
tetha0=0.75
tetha1=1.2
tetha2=1.5
tetha3=0.5
tetha4=1.25
tetha5=1.55
tetha6=0.25
tetha7=1.0
y0 = [1.0,1.0]
y0polet=[0,0,V0,tetha0]
y0polet1=[0,0,V0,tetha1]
y0polet2=[0,0,V0,tetha2]
y0polet3=[0,0,V0,tetha3]
y0polet4=[0,0,V0,tetha4]
y0polet5=[0,0,V0,tetha5]
y0polet6=[0,0,V0,tetha6]
y0polet7=[0,0,V0,tetha7]
def rungekutta2(f, y0polet, t, args=()):
    n = len(t)
    y = np.zeros((n, len(y0polet)))
    y[0] = y0polet
    for i in range(n - 1):
        h = t[i+1] - t[i]
        y[i+1] = y[i] + h * f(y[i] + f(y[i], t[i], *args) * h / 2., t[i] + h / 2., *args)
    return y
'''
err1=[]
err2=[]
sol4=[]
ytoch=[]
j=0
toch=[]
for i in range(1,978,10):
    toch.append(i)
    t4 = np.linspace(0, 4, i)
    sol4 = rungekutta2(test, y0, t4)
    xtoch = np.linspace(0, 4, i)
    ytoch = [[exp(j),exp(-j)] for j in xtoch]
    num=[abs(sol4[j][0]-ytoch[j][0])]
    err1.append(max(num))
    err2.append(max(num)/pow(i,4))
    j+=1
plt.plot(toch[::-1], err1, label='решение е')
plt.plot(toch[::-1], err2, label='решение e/h^4')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
'''
'''
t4 = np.linspace(0, 4, 978)
sol4 = rungekutta2(test, y0, t4)
xtoch = np.linspace(0, 4, 10)
ytoch1 = [[exp(j)] for j in xtoch]
ytoch2 = [[exp(-j)] for j in xtoch]
plt.plot(t4, sol4[:, 0], label='численное решение - e^x')
plt.plot(t4, sol4[:, 1], label='численное решение - e^(-x)')
plt.plot(xtoch,ytoch1[0:],label='точное решение - e^x')
plt.plot(xtoch,ytoch2[0:],label='точное решение - e^(-x)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
'''
t=np.linspace(0,7.67,100)
sol1=rungekutta2(polet, y0polet, t, args=(m, C, rho, S, g, V0, tetha0))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=0.75$')

t=np.linspace(0,10.35,100)
sol1=rungekutta2(polet, y0polet1, t, args=(m, C, rho, S, g, V0, tetha1))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=1.2$')

t=np.linspace(0,11.05,100)
sol1=rungekutta2(polet, y0polet2, t, args=(m, C, rho, S, g, V0, tetha2))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=1.5$')
t=np.linspace(0,5.5,100)
sol1=rungekutta2(polet, y0polet3, t, args=(m, C, rho, S, g, V0, tetha3))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=0.5$')
t=np.linspace(0,10.525,100)
sol1=rungekutta2(polet, y0polet4, t, args=(m, C, rho, S, g, V0, tetha4))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=1.25$')
t=np.linspace(0,11.07,100)
sol1=rungekutta2(polet, y0polet5, t, args=(m, C, rho, S, g, V0, tetha5))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=1.55$')
t=np.linspace(0,2.9,100)
sol1=rungekutta2(polet, y0polet6, t, args=(m, C, rho, S, g, V0, tetha6))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=0.25$')
t=np.linspace(0,9.35,100)
sol1=rungekutta2(polet, y0polet7, t, args=(m, C, rho, S, g, V0, tetha7))
plt.plot(sol1[:,0], sol1[:,1], label=r'$\theta_0=1.0$')
plt.legend(loc='best')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.show()

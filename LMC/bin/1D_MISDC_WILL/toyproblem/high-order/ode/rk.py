from math import *

quad_pts = [0, (1-sqrt(1/5.0))*0.5, (1+sqrt(1/5.0))*0.5, 1]

def interp(f, dt, t):
    return (2*pow(dt,3)*f[0] + 10*pow(t,3)*(-f[0] + sqrt(5)*f[1] - sqrt(5)*f[2] + f[3]) + 
    pow(dt,2)*t*(-12*f[0] + 5*(1 + sqrt(5))*f[1] + 5*f[2] - 5*sqrt(5)*f[2] + 2*f[3]) - 
    5*dt*pow(t,2)*(-4*f[0] + f[1] + 3*sqrt(5)*f[1] + f[2] - 3*sqrt(5)*f[2] + 2*f[3]))/(2.*pow(dt,3))


def rk4(f, x0, y0, x1, n):
    vx = [0]*(n + 1)
    vy = [0]*(n + 1)
    h = (x1 - x0)/float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        print x
        k1 = h*f(x, y)
        k2 = h*f(x + 0.5*h, y + 0.5*k1)
        k3 = h*f(x + 0.5*h, y + 0.5*k2)
        k4 = h*f(x + h, y + k3)
        vx[i] = x = x0 + i*h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4)/6.0
    return vx, vy
 
def f(x, y):
    return interp(quad_pts, 1, x) + y
    #return x + y
 
vx, vy = rk4(f, 0, 0, 1, 10)

print 'rk: ', vy[-1]
print 'ex: ', exp(1)-2

#for x, y in list(zip(vx, vy))[::10]:
#    print(x, y, y - (4 + x*x)**2/16)

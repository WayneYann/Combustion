import pde3
import finalpde
import matplotlib.pyplot as plt

#a = -1.0
#eps = 0.4
#r = -7.0

#a = -1.0
#eps = 30.0
#r = -70.0

a = -1.0
eps = 10.0
r = -300.0

Nx = 300
dt = (20.0/Nx)/5.0
j = 0
#FT = dt*200
FT = 2.0

(x1, y1) = pde3.solve_pde(a, eps, r, Nx=Nx*2**(j), dt=dt/(2**(2*j)),
                        FT=FT, max_iter=2)
(x, y) = finalpde.solve_pde(a, eps, r, Nx=Nx*2**(j), dt=dt/(2**(2*j)),
                        FT=FT, max_iter=2)

plt.cla()
plt.plot(x,y, x1, y1)
plt.ylim((-0.1,1.1))
plt.show(block=True)

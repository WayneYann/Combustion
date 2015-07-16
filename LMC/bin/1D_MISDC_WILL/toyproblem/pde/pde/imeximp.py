#    # simple IMEX
#    if IMEX:
#        y[n+1][1:-1] = np.linalg.solve(I - dt*eps*Dxx,
#                                         y[n][1:-1]
#                                       + dt*(np.dot(A, y[n][1:-1])
#                                       + (-a/(2*h) + eps/h**2)*bc))
#        return
#    # purely implicit method
#    elif IMP:
#        y[n+1][1:-1] = np.linalg.solve(I - dt*(a*Dx + eps*Dxx),
#                                         y[n][1:-1]
#                                       + dt*(-a/(2*h) + eps/h**2)*bc)
#        return

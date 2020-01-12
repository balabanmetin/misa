
from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint, Bounds, BFGS
import numpy as np
from misa.Method import FM, OLS

MIN_X = 1e-5
MAX_X = 5.0
maxIter=5000



def optimize_for_two(branch1, branch2, tree, obs_dist, model_name, method_name):


    dlist=[tree.distance_between(branch1,branch2)]
    if branch1.parent!= tree.root:
        dlist += [tree.distance_between(branch1.parent, branch2)]
    if branch2.parent != tree.root:
        dlist += [tree.distance_between(branch1, branch2.parent)]
    d = max(dlist)
    d=0.0219976


    mvec = [obs_dist[k] for k in sorted(obs_dist)]
    n=len(mvec)
    x0=np.array(mvec+mvec+[MIN_X,MIN_X,MIN_X,MIN_X])

    if model_name == "HAR":
        # harmonic as y approaches to infinity , harmonic mean of x and y is 2*x, thus the lower bound
        bounds = Bounds(np.array(mvec+mvec+[0,0,0,0])/2,
                        np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))


        def cons_f(x):
            return [2*x[k]*x[n+k] -mvec[k]*x[k] - mvec[k]*x[k+n] for k in range(n)]

        def cons_g(x):
            res = [[0]*(2*n+4) for k in range(n)]
            for k in range(n):
                res[k][k] = 2*x[n+k] - mvec[k]
                res[k][k+n] = 2*x[k] - mvec[k]
            return res

        constraint = NonlinearConstraint(cons_f, 0, 0, jac=cons_g, hess='3-point')

    elif model_name == "JAC":

        hid=np.array(mvec+mvec)
        lb = 1 - (1.5)**(1/31)*(1-hid)

        bounds = Bounds(np.concatenate((lb, np.array([0]*4)), axis=None),
                        np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))


        def cons_f(x):
            x1, x2, x3, x4 = x[-4:]
            return [2*((1-x[k])**31 + (1-x[k+n])**31 - (1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**31)
                    - (1-mvec[k])**31*(3-(1-(d-x1-x3+x2+x4))**31) for k in range(n)]

        def cons_g(x):
            x1, x2, x3, x4 = x[-4:]
            res = [[0]*(2*n+4) for k in range(n)]
            for k in range(n):
                res[k][k]   = -2*31*(1-x[k])**30 + 31/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**30
                res[k][k+n] = -2*31*(1-x[k+n])**30 + 31/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**30
                temp = - 31/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**30 + 31*(1-mvec[k])**31*(1-(d-x1-x3+x2+x4))**30
                res[k][-4] = temp
                res[k][-2] = temp
                res[k][-3] = -temp
                res[k][-1] = -temp
            return res

        constraint = NonlinearConstraint(cons_f, 0, 0, jac=cons_g, hess='cs')

    elif model_name == "SIMPJAC":

        hid=np.array(mvec+mvec)
        lb = 1 - (1.5)**(1/31)*(1-hid)
        lb_nonnegative = np.array([ max(i,0) for i in lb])
        bounds = Bounds(np.concatenate((lb_nonnegative, np.array([0]*4)), axis=None),
                        np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))


        def cons_f(x):
            return [2*((1-x[k])**31 + (1-x[k+n])**31 - (1 - (x[k]+x[k+n]+d)/2)**31)
                    - (1-mvec[k])**31*(3-(1-d)**31) for k in range(n)]

        def cons_g(x):
            res = [[0]*(2*n+4) for k in range(n)]
            for k in range(n):
                res[k][k]   = -2*31*(1-x[k])**30 + 2*31/2*(1 - (x[k]+x[k+n]+d)/2)**30
                res[k][k+n] = -2*31*(1-x[k+n])**30 + 2*31/2*(1 - (x[k]+x[k+n]+d)/2)**30
            return res

        def cons_h_v(x,v):
            res = [[0] * (2 * n + 4) for k in range(2*n+4)]
            for k in range(n):
                res[k][k] = v[k] * (2*31*30*(1-x[k])**29 - 2*31/2*30/2*(1 - (x[k]+x[k+n]+d)/2)**29)
                res[k+n][k+n] = v[k] * (2*31*30*(1-x[k+n])**29 - 2*31/2*30/2*(1 - (x[k]+x[k+n]+d)/2)**29)
                res[k][k + n] = res[k+n][k] = v[k] * (- 2*31/2*30/2*(1 - (x[k]+x[k+n]+d)/2)**29)
            return res

        constraint = NonlinearConstraint(cons_f, 0, 0, jac=cons_g, hess=cons_h_v)


        A=[[0.0] * (2 * n + 4) for i in range(3*n)]
        ub = [0.0]*(3*n)
        lb = [0.0]*(3*n)
        for i in range(n):
            A[i][i]=A[i][i+n]=1
            lb[i]=d
            ub[i]=MAX_X

            A[i+n][i] = 1
            A[i + n][i+n] = -1
            lb[i+n]=-MAX_X
            ub[i+n]=d

            A[i + 2*n][i] = -1
            A[i + 2*n][i + n] = 1
            lb[i + 2*n] = -MAX_X
            ub[i + 2*n] = d

        constraint_triangle = LinearConstraint(A,lb=lb,ub=ub)
    else: #linear model
        bounds = Bounds(np.array([0]*(2*n+4)),np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))

        # constraints depend on model too
        cons_mtr = [[0.0] * (2 * n + 4) for i in range(n)]
        for i in range(n):
            cons_mtr[i][i] = 0.5
            cons_mtr[i][i + n] = 0.5
        constraint = LinearConstraint(cons_mtr, mvec, mvec)


    if method_name == "FM":
        f = FM.f
        g = FM.g
        h = FM.h
        h_p = FM.h_p
        result = minimize(fun=f, method="trust-constr", x0=x0, bounds=bounds, args=(branch1, branch2),
                          constraints=[constraint,constraint_triangle],
                          options={'disp': True, 'verbose': 1, 'maxiter': maxIter}, jac=g, hessp='3-point')
    else:
        f = OLS.f
        g = OLS.g
        h = OLS.h
        h_p = OLS.h_p

        try:
            result = minimize(fun=f, method="trust-constr", x0=x0, bounds=bounds, args=(branch1, branch2), constraints=[constraint,constraint_triangle],
                      options={'disp': True, 'verbose': 1, 'maxiter': maxIter} , jac=g, hessp=h_p )
        except Exception as e:
            return (None, branch1, branch2)

    return (result, branch1, branch2)

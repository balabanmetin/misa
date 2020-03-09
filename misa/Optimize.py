
from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint, Bounds, BFGS
import numpy as np
from misa.Method import FM, OLS

MIN_X = 0
MAX_X = 0.67



def optimize_for_two(branch1, branch2, tree, obs_dist, model_name, method_name,maxIter,ks):


    dlist=[tree.distance_between(branch1,branch2)]
    if branch1.parent!= tree.root:
        dlist += [tree.distance_between(branch1.parent, branch2)]
    if branch2.parent != tree.root:
        dlist += [tree.distance_between(branch1, branch2.parent)]
    d = max(dlist)

    mvec = [obs_dist[k] for k in sorted(obs_dist)]
    n=len(mvec)
    obs_min = 1-(2/(3-(1-d)**ks))**(1.0/ks)
    for i in range(n):
        mvec[i]=max(mvec[i],obs_min+0.001)

    def init_zero_obj():
        x0=np.array(mvec+mvec+[MIN_X,MIN_X,MIN_X,MIN_X])
        for k, v in branch1.Rd.items():
            x0[k] = v + branch1.edge_length
        for k, v in branch1.Sd.items():
            x0[k] = v
        for k, v in branch2.Rd.items():
            x0[k + n] = v + branch2.edge_length
        for k, v in branch2.Sd.items():
            x0[k + n] = v
            return x0

    def init_zero_constr_viol():
        nmvec = 1- (1-np.array(mvec))*((3-(1-d)**ks)/2)**(1/ks)
        hypo1=[0]*n
        hypo2=[0]*n
        res = [0]*(2*n+4)
        for k, v in branch1.Rd.items():
            hypo1[k] = v + branch1.edge_length
        for k, v in branch1.Sd.items():
            hypo1[k] = v
        for k, v in branch2.Rd.items():
            hypo2[k] = v + branch2.edge_length
        for k, v in branch2.Sd.items():
            hypo2[k] = v
        for i in range(n):
            if hypo1[i] <= hypo2[i]:
                res[i] = nmvec[i]
                res[i+n] = nmvec[i]+d
            else:
                res[i+n] = nmvec[i]
                res[i] = nmvec[i] + d
        return res

    def init_lazy():
        return np.array(mvec+mvec+[MIN_X,MIN_X,MIN_X,MIN_X])

    #x0=init_zero_obj()
    x0=init_zero_constr_viol()
    #x0=init_lazy()

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
        lb = 1 - (1.5)**(1/ks)*(1-hid)

        bounds = Bounds(np.concatenate((lb, np.array([0]*4)), axis=None),
                        np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))


        def cons_f(x):
            x1, x2, x3, x4 = x[-4:]
            return [2*((1-x[k])**ks + (1-x[k+n])**ks - (1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**ks)
                    - (1-mvec[k])**ks*(3-(1-(d-x1-x3+x2+x4))**ks) for k in range(n)]

        def cons_g(x):
            x1, x2, x3, x4 = x[-4:]
            res = [[0]*(2*n+4) for k in range(n)]
            for k in range(n):
                res[k][k]   = -2*ks*(1-x[k])**(ks-1) + ks/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**(ks-1)
                res[k][k+n] = -2*ks*(1-x[k+n])**(ks-1) + ks/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**(ks-1)
                temp = - ks/2*(1 - (x[k]+x[k+n]+d-x1-x3+x2+x4)/2)**(ks-1) + ks*(1-mvec[k])**ks*(1-(d-x1-x3+x2+x4))**(ks-1)
                res[k][-4] = temp
                res[k][-2] = temp
                res[k][-3] = -temp
                res[k][-1] = -temp
            return res

        constraint = NonlinearConstraint(cons_f, 0, 0, jac=cons_g, hess='cs')

    elif model_name == "SIMPJAC":

        hid=np.array(mvec+mvec)
        lb = 1 - (1.5)**(1/ks)*(1-hid)
        lb_nonnegative = np.array([ max(i,0) for i in lb])
        bounds = Bounds(np.concatenate((lb_nonnegative, np.array([0]*4)), axis=None),
                        np.array([MAX_X]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))


        def cons_f(x):
            return [2*((1-x[k])**ks + (1-x[k+n])**ks - (1 - (x[k]+x[k+n]+d)/2)**ks)
                    - (1-mvec[k])**ks*(3-(1-d)**ks) for k in range(n)]

        def cons_g(x):
            res = [[0]*(2*n+4) for k in range(n)]
            for k in range(n):
                res[k][k]   = -2*ks*(1-x[k])**(ks-1) + 2*ks/2*(1 - (x[k]+x[k+n]+d)/2)**(ks-1)
                res[k][k+n] = -2*ks*(1-x[k+n])**(ks-1) + 2*ks/2*(1 - (x[k]+x[k+n]+d)/2)**(ks-1)
            return res

        def cons_h_v(x,v):
            res = [[0] * (2 * n + 4) for k in range(2*n+4)]
            for k in range(n):
                res[k][k] = v[k] * (2*ks*(ks-1)*(1-x[k])**(ks-2) - 2*ks/2*(ks-1)/2*(1 - (x[k]+x[k+n]+d)/2)**(ks-2))
                res[k+n][k+n] = v[k] * (2*ks*(ks-1)*(1-x[k+n])**(ks-2) - 2*ks/2*(ks-1)/2*(1 - (x[k]+x[k+n]+d)/2)**(ks-2))
                res[k][k + n] = res[k+n][k] = v[k] * (- 2*ks/2*(ks-1)/2*(1 - (x[k]+x[k+n]+d)/2)**(ks-2))
            return res

        constraint = NonlinearConstraint(cons_f, 0, 0, jac=cons_g, hess=cons_h_v)
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
                          constraints=[constraint],
                          options={'disp': True, 'verbose': 1, 'maxiter': maxIter}, jac=g, hessp='3-point')
    else:
        f = OLS.f
        g = OLS.g
        h = OLS.h
        h_p = OLS.h_p

        try:
            result = minimize(fun=f, method="trust-constr", x0=x0, bounds=bounds, args=(branch1, branch2), constraints=[constraint],
                      options={'disp': True, 'verbose': 1, 'maxiter': maxIter} , jac=g, hess=h )
        except Exception as e:
            return (None, branch1, branch2)
    #print(result.fun , branch1.edge_index, branch2.edge_index)
    return (result, branch1, branch2)


from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint, Bounds, BFGS
import numpy as np
from misa.Method import FM, OLS
import math

MIN_X = 0
MAX_X = 0.67



def optimize_for_two(branch1, branch2, tree, obs_containment, model_name, method_name, maxIter, k, alpha):


    dlist=[tree.distance_between(branch1,branch2)]
    if branch1.parent!= tree.root:
        dlist += [tree.distance_between(branch1.parent, branch2)]
    if branch2.parent != tree.root:
        dlist += [tree.distance_between(branch1, branch2.parent)]
    d = max(dlist)

    mvec = [obs_containment[k] for k in sorted(obs_containment)]
    n=len(mvec)

    def init_zero_obj():
        x0=np.array(mvec+mvec+[MIN_X,MIN_X,MIN_X,MIN_X])
        for key, v in branch1.Rd.items():
            x0[key] = v + branch1.edge_length
        for key, v in branch1.Sd.items():
            x0[key] = v
        for key, v in branch2.Rd.items():
            x0[key + n] = v + branch2.edge_length
        for key, v in branch2.Sd.items():
            x0[key + n] = v
        for i in range(n):
            x0[i] = math.exp(-k*x0[i])
            x0[i+n] = math.exp(-k*x0[i+n])

        return x0

    def init_zero_constr_viol():
        nmvec = 1- (1-np.array(mvec))*((3-(1-d)**31)/2)**(1/31)
        hypo1=[0]*n
        hypo2=[0]*n
        res = np.array([0]*(2*n+4))
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

    x0=init_zero_obj()
    #x0=init_zero_constr_viol()
    #x0=init_lazy()


    lb = np.array([math.exp(-k*MAX_X)]*(2*n))
    lb=np.concatenate((lb, np.array([0] * 4)), axis=None)
    ub=np.array([1]*(2*n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X])
    bounds = Bounds(lb,ub)


    if method_name == "FM":
        f = FM.f
        g = FM.g
        h = FM.h
        h_p = FM.h_p
        result = minimize(fun=f, method="trust-constr", x0=x0, bounds=bounds, args=(branch1, branch2),
                          options={'disp': True, 'verbose': 1, 'maxiter': maxIter}, jac=g, hessp='3-point')
    else:
        f = OLS.f
        g = OLS.g
        h = OLS.h
        h_p = OLS.h_p

        try:
            result = minimize(fun=f, method="trust-constr", x0=x0,  args=(branch1, branch2, alpha, k, d, mvec),
                      options={'disp': True, 'verbose': 1, 'maxiter': maxIter} , jac=g,
                              constraints=[LinearConstraint(A=np.eye(2*n+4), lb=lb, ub=ub)], hess='3-point' )
        except Exception as e:
            print("optimization failed")
            import pdb; pdb.set_trace()
            return (None, branch1, branch2)
    print(result.fun , branch1.edge_index, branch2.edge_index)
    for i in range(2*n):
        print(-1/k*math.log(result.x[i]),)
    print("")
    return (result, branch1, branch2)

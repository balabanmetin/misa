
from abc import ABC ,abstractmethod
import numpy as np

# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with \delta^2 scaling
# BE: Beyer --Least Squares with \delta scaling


class Method():
    @staticmethod
    def f(x, *args):
        pass

    @staticmethod
    def g(x, *args):
        pass

    @staticmethod
    def h(x, *args):
        pass


class OLS(Method):

    @staticmethod
    def f(x, *args):
        branch1 = args[0]
        branch2 = args[1]
        w = args[2]
        signs=args[3]
        xAB, x1, x2, x3, x4 = x[-5:]
        acc = 0
        n = int((len(x) - 5) / 2)
        for k, v in branch1.Rd.items():
            acc += (x[k] - (v + branch1.edge_length - x1 + x2)) ** 2
        for k, v in branch1.Sd.items():
            acc += (x[k] - (v + x1 + x2)) ** 2

        for k, v in branch2.Rd.items():
            acc += (x[k + n] - (v + branch2.edge_length - x3 + x4)) ** 2
        for k, v in branch2.Sd.items():
            acc += (x[k + n] - (v + x3 + x4)) ** 2

        # w is distance between distal tips of branch1 and branch2
        #(xAB - w + x1 - x2 + x3 - x4)
        acc += (-w  + np.dot(np.array([xAB,x1,x2,x3,x4]),signs))**2
        return acc

    @staticmethod
    def g(x, *args):
        branch1 = args[0]
        branch2 = args[1]
        w = args[2]
        signs=args[3]
        x1, x2, x3, x4 = x[-4:]
        xAB = x[-5]
        res = len(x) * [0]
        n = int((len(x) - 5) / 2)

        for k, v in branch1.Rd.items():
            res[k] += 2 * (x[k] - (v + branch1.edge_length - x1 + x2))
            res[-4] += 2 * (x[k] - (v + branch1.edge_length - x1 + x2))
            res[-3] += - 2 * (x[k] - (v + branch1.edge_length - x1 + x2))

        for k, v in branch1.Sd.items():
            res[k] += 2 * (x[k] - (v + x1 + x2))
            res[-4] += - 2 * (x[k] - (v + x1 + x2))
            res[-3] += - 2 * (x[k] - (v + x1 + x2))

        for k, v in branch2.Rd.items():
            res[k + n] += 2 * (x[k + n] - (v + branch2.edge_length - x3 + x4))
            res[-2] += 2 * (x[k + n] - (v + branch2.edge_length - x3 + x4))
            res[-1] += - 2 * (x[k + n] - (v + branch2.edge_length - x3 + x4))

        for k, v in branch2.Sd.items():
            res[k + n] += 2 * (x[k + n] - (v + x3 + x4))
            res[-2] += - 2 * (x[k + n] - (v + x3 + x4))
            res[-1] += - 2 * (x[k + n] - (v + x3 + x4))

        #res[-5] = 2 * (xAB - (np.dot(np.array([w,x1,x2,x3,x4]),signs)))
        res[-5:] += 2*(-w  + np.dot(np.array([xAB,x1,x2,x3,x4]),signs))*signs

        # res[-4] += -signs[-4]* res[-5]
        # res[-3] += -signs[-3]* res[-5]
        # res[-2] += -signs[-2]* res[-5]
        # res[-1] += -signs[-1]* res[-5]
        return res

    @staticmethod
    def h(x, *args):
        branch1 = args[0]
        branch2 = args[1]
        signs = args[3]
        n = int((len(x) - 5) / 2)
        H = np.diag([2] * len(x), 0)

        H[-4][-4] = H[-3][-3] = H[-2][-2] = H[-1][-1]= 2 * n
        H[-4][-3] = H[-3][-4] = -len(branch1.Rd)*2 + len(branch1.Sd)*2
        H[-1][-2] = H[-2][-1] = -len(branch2.Rd)*2 + len(branch2.Sd)*2

        for k, v in branch1.Rd.items():
            H[k][-4] = H[-4][k] = 2

        for k, v in branch1.Sd.items():
            H[k][-4] = H[-4][k] = -2

        for k, v in branch2.Rd.items():
            H[k+n][-2] = H[-2][k+n] = 2

        for k, v in branch2.Sd.items():
            H[k+n][-2] = H[-2][k+n] = -2

        H[0:n, -3] = H[n:2*n, -1] = H[-1, n:2*n] = H[-3, 0:n] = -2


        H[-5][-5]=0
        last_term_hessian = 2*np.outer(signs,signs)
        for i in range(5):
            for j in range(5):
                H[2*n+i][2*n+j] += last_term_hessian[i][j]
        return H

    @staticmethod
    def h_p(x, p, *args):
        n = int((len(x) - 4) / 2)
        res = np.zeros_like(x)
        res[0:n]  = 4*p[0:n]   - 4*p[-1]
        res[n:2*n]= 4*p[n:2*n] - 4*p[-3]
        res[-4:]  = 4*p[-4:]
        res[-3] += -4*np.sum(p[n:2*n])
        res[-1] += -4*np.sum(p[0:n])
        return res

class FM(Method):

    @staticmethod
    def f(x, *args):
        branch1 = args[0]
        branch2 = args[1]
        x1, x2, x3, x4 = x[-4:]
        acc = 0
        n = int((len(x) - 4) / 2)
        for k, v in branch1.Rd.items():
            acc += (1 - (v + branch1.edge_length - x1 + x2)/x[k]) ** 2
        for k, v in branch1.Sd.items():
            acc += (1 - (v + x1 + x2)/x[k]) ** 2

        for k, v in branch2.Rd.items():
            acc += (1 - (v + branch2.edge_length - x3 + x4)/x[k]) ** 2
        for k, v in branch2.Sd.items():
            acc += (1 - (v + x3 + x4)/x[k]) ** 2

        return acc

    @staticmethod
    def g(x, *args):
        branch1 = args[0]
        branch2 = args[1]
        x1, x2, x3, x4 = x[-4:]
        res = len(x) * [0]
        n = int((len(x) - 4) / 2)

        for k, v in branch1.Rd.items():
            res[k] += 2 * (v + branch1.edge_length - x1 + x2) * (x[k] - (v + branch1.edge_length - x1 + x2))/x[k]**3
            res[-4] += 2 * (x[k] - (v + branch1.edge_length - x1 + x2))/x[k]**2
            res[-3] += - 2 * (x[k] - (v + branch1.edge_length - x1 + x2))/x[k]**2

        for k, v in branch1.Sd.items():
            res[k] += 2 * (v + x1 + x2)*(x[k] - (v + x1 + x2))/x[k]**3
            res[-4] += - 2 * (x[k] - (v + x1 + x2))/x[k]**2
            res[-3] += - 2 * (x[k] - (v + x1 + x2))/x[k]**2

        for k, v in branch2.Rd.items():
            res[k + n] += 2 * (v + branch2.edge_length - x3 + x4) * (x[k + n] - (v + branch2.edge_length - x3 + x4))/x[k + n]**3
            res[-2] += 2 * (x[k + n] - (v + branch2.edge_length - x3 + x4))/x[k+n]**2
            res[-1] += - 2 * (x[k + n] - (v + branch2.edge_length - x3 + x4))/x[k+n]**2

        for k, v in branch2.Sd.items():
            res[k + n] += 2 * (v + x3 + x4)* (x[k + n] - (v + x3 + x4))/x[k + n]**3
            res[-2] += - 2 * (x[k + n] - (v + x3 + x4))/x[k+n]**2
            res[-1] += - 2 * (x[k + n] - (v + x3 + x4))/x[k+n]**2

        return res

    @staticmethod
    def h(x, *args):

        H = np.diag([4] * len(x), 0)
        H[:, -1] = H[:, -3] = H[-3, :] = H[-1, :] = 4

        return H

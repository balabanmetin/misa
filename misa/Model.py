
import numpy as np
from abc import ABC, abstractmethod


class Model(ABC):

    @abstractmethod
    def bounds_and_const(self, mvec, branch1, branch2):
        pass



class Linear_model(Model):
    def bounds_and_const(self, mvec, branch1, branch2):
        bounds = Bounds(np.array([0] * (2 * n + 4)),
                        np.array([MAX_X] * (2 * n) + [branch1.edge_length, MAX_X, branch2.edge_length, MAX_X]))

        # constraints depend on model too
        cons_mtr = [[0.0] * (2 * n + 4) for i in range(n)]
        for i in range(n):
            cons_mtr[i][i] = 0.5
            cons_mtr[i][i + n] = 0.5
        linear_constraint = LinearConstraint(cons_mtr, mvec, mvec)
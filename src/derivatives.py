"""
Different classes for finite difference derivatives of a function with several variables
"""

import numpy as np

class FirstOrder:
    def __init__(self):
        pass
    
    def differentiateFirstFO(self, x, y, h):
        """
        Compute first order first derivative using finite difference method
        """

        data_dim = x.shape[-1]
        diff = np.empty(data_dim)

        diff[0]  = (y[1] - y[0])/h
        diff[-1] = (y[-1] - y[-2])/h

        for i in range(1,data_dim-1):
            diff[i] = (y[i+1] - y[i-1])/(2*h)

        return diff

    def differentiateSecondFO(self, x):
        """
        Compute first order second derivative using finite difference method
        """

    def differentiateThirdFO(self, x):
        """
        Compute first order third derivative using finite difference method
        """

    def differentiateFourthFO(self, x):
        """
        Compute first order fourth derivative using finite difference method
        """

class SecondOrder:
    def __init__(self):
        pass
    
    def differentiateCentralFirstSO(self, x):
        """
        Compute second order central first derivative using finite difference method
        """

    def differentiateCentralSecondSO(self, x):
        """
        Compute second order central second derivative using finite difference method
        """

    def differentiateCentralThirdSO(self, x):
        """
        Compute second order central third derivative using finite difference method
        """

    def differentiateCentralFourthSO(self, x):
        """
        Compute second order central fourth derivative using finite difference method
        """

    def differentiateForwardFirstSO(self, x):
        """
        Compute second order forward first derivative using finite difference method
        """

    def differentiateForwardSecondSO(self, x):
        """
        Compute second order forward second derivative using finite difference method
        """

    def differentiateForwardThirdSO(self, x):
        """
        Compute second order forward third derivative using finite difference method
        """

    def differentiateForwardFourthSO(self, x):
        """
        Compute second order forward fourth derivative using finite difference method
        """

    def differentiateBackwardFirstSO(self, x):
        """
        Compute second order backward first derivative using finite difference method
        """

    def differentiateBackwardSecondSO(self, x):
        """
        Compute second order backward second derivative using finite difference method
        """

    def differentiateBackwardThirdSO(self, x):
        """
        Compute second order backward third derivative using finite difference method
        """

    def differentiateBackwardFourthSO(self, x):
        """
        Compute second order backward fourth derivative using finite difference method
        """
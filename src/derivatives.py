"""
Different classes for finite difference derivatives of a function with several variables
"""

import numpy as np

class FirstOrder:
    def __init__(self):
        pass
    
    def differentiateFirstFO(self, x, y, h):
        """
        Compute first order first derivative using finite difference method.
        """
        num_variables = x.shape[0]
        data_dim      = x.shape[-1]
        diff          = np.empty([num_variables, data_dim])

        for j in range(num_variables):

            diff[j, 0]  = (y[j,1] - y[j,0])/h
            diff[j, -1] = (y[j,-1] - y[j,-2])/h

            for i in range(1,data_dim-1):
                diff[j,i] = (y[j,i+1] - y[j, i-1])/(2*h)

        return diff

    def differentiateSecondFO(self, x, y, h):
        """
        Compute first order second derivative using finite difference method
        """

        num_variables = x.shape[0]
        data_dim      = x.shape[-1]
        diff          = np.empty([num_variables, data_dim])

        for j in range(num_variables):

            diff[j, 0]  = (y[j,2] - 2*y[j, 1] + y[j,0])/h**2
            diff[j, -1] = (y[j,-1] - 2*y[j,-2] + y[j, -3])/h**2

            for i in range(1,data_dim-1):
                diff[j,i] = (y[j,i+1] - 2*y[j, i] +  y[j, i-1])/h**2

        return diff

    def differentiateThirdFO(self, x, y, h):
        """
        Compute first order third derivative using finite difference method
        """

        num_variables = x.shape[0]
        data_dim      = x.shape[-1]
        diff          = np.empty([num_variables, data_dim])

        for j in range(num_variables):

            diff[j, 0]  = (y[j,3] - 3*y[j, 2] + 3*y[j, 1] - y[j,0])/h**3
            diff[j, 1]  = (y[j,4] - 3*y[j, 3] + 3*y[j, 2] - y[j,1])/h**3

            diff[j, -1] = (y[j,-1] - 3*y[j,-2] + 3*y[j, -3] - y[j, -4])/h**3
            diff[j, -2] = (y[j,-2] - 3*y[j,-3] + 3*y[j, -4] - y[j, -5])/h**3

            for i in range(2,data_dim-2):
                diff[j,i] = (y[j,i+2] - 2*y[j,i+1] + 2*y[j, i-1] -y[j,i-2])/(2*(h**3))

        return diff


    def differentiateFourthFO(self, x, y, h):
        """
        Compute first order fourth derivative using finite difference method
        """

        num_variables = x.shape[0]
        data_dim      = x.shape[-1]
        diff          = np.empty([num_variables, data_dim])

        for j in range(num_variables):

            diff[j, 0]  = (y[j,4] - 4*y[j, 3] + 6*y[j, 2] - 4*y[j,1] + y[j,0])/h**4
            diff[j, 1]  = (y[j,5] - 4*y[j, 4] + 6*y[j, 3] - 4*y[j,2] + y[j,1])/h**4
            diff[j, 2]  = (y[j,6] - 4*y[j, 5] + 6*y[j, 4] - 4*y[j,3] + y[j,2])/h**4

            diff[j, -1] = (y[j,-1] - 4*y[j,-2] + 6*y[j, -3] - 4*y[j,-4] + y[j,-5])/h**4
            diff[j, -2] = (y[j,-2] - 4*y[j,-3] + 6*y[j, -4] - 4*y[j,-5] + y[j,-6])/h**4
            diff[j, -3] = (y[j,-3] - 4*y[j,-4] + 6*y[j, -5] - 4*y[j,-6] + y[j,-7])/h**4

            for i in range(3,data_dim-2):
                diff[j,i] = (y[j,i+2] - 4*y[j,i+1] + 6*y[j,i] - 4*y[j, i-1] +  y[j, i-2])/h**4

        return diff

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
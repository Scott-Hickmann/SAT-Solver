from decimal import Decimal, getcontext
import numpy as np

class ODESolver:
    def __init__(self, func, precision=4):
        """
        Initializes the ODE solver.
        
        Parameters:
        - func: a function of the form f(t, y) representing the ODE dy/dt = f(t, y).
        """
        self.func = func
        getcontext().prec = precision

    def euler_step(self, t, y, dt):
        """
        Perform one step of the Euler method to solve the ODE.
        
        Parameters:
        - t: current time
        - y: current value
        - dt: time step size
        
        Returns:
        - new value of y after time dt
        """
        return y + dt * self.func(t, y)

    def solve(self, y0, t0, t_end, dt):
        """
        Solve the ODE using the Euler method from time t0 to t_end with initial value y0.
        
        Parameters:
        - y0: initial value
        - t0: start time
        - t_end: end time
        - dt: time step size
        
        Returns:
        - times: list of time points
        - values: list of y values corresponding to the time points
        """
        times = [t0]
        values = [y0]
        
        t = t0
        y = y0
        
        while t < t_end:
            y = self.euler_step(t, y, dt)
            t += dt
            
            times.append(t)
            values.append(y)
        
        return np.array(times), np.array(values)

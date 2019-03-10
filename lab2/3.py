import matplotlib.pyplot as plt
import numpy as np
from math import pi

def P(t, k):
    return np.arctan(t/4-1/2) + pi*k

def SIM(P, Iterations, k):
                                                                                                                                                       
    x = 1
                                                                                                                                                                               
    for i in range(Iterations):        
        x = P(x, k)        
                
    return x

k = 0

print(f'k-й корень примерно равен: {round(SIM(P, 10, k), 4)}')

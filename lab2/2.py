import matplotlib.pyplot as plt
import numpy as np
from math import pi

plt.figure(1)

plt.subplot(111)

x = np.arange(-20, 20, 0.01)
y = [-1.00, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1.00]

def df(t):
    return 4/(16 + (t-2)**2)

plt.plot(x, df(x), 'b-', x, np.zeros(len(x)), 'k-', np.zeros(len(x)), 10*x, 'k-')


plt.xlim(-20, 20)
plt.ylim(-1, 1)
plt.grid(True)
plt.show()

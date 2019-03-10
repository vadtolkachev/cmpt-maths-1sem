import matplotlib.pyplot as plt
import numpy as np
from math import pi

plt.figure(1)

plt.subplot(121)


def f1(t):
	return np.tan(t)

def f2(t):
	return t/4 - 1/2

for i in range(-3,4):
	x = np.arange(-pi/2 + 0.01 + i*pi, pi/2 - 0.01 + i*pi, 0.01)
	plt.plot(x, f1(x), 'b-', x, np.zeros(len(x)), 'k-', np.zeros(len(x)), 10*x, 'k-')

x = np.arange(-10, 10, 0.01)
plt.plot(x, f2(x), 'b-', x, np.zeros(len(x)), 'k-', np.zeros(len(x)), 10*x, 'k-')

plt.annotate('y', xy=(0, 5), xytext=(0, 4),
            arrowprops=dict(facecolor='black', shrink=0.5),
            )
plt.annotate('', xy=(10.0, 0), xytext=(9.9, 0),
            arrowprops=dict(facecolor='black', shrink=0.5),
            )
plt.text(8.5, 0.5, r'x')

plt.xlim(-10,10)
plt.ylim(-5,5)
plt.grid(True)

plt.subplot(122)
plt.plot(x, f1(x), 'b-', x, np.zeros(len(x)), 'k-', np.zeros(len(x)), 10*x, 'k-')
plt.plot(x, f2(x), 'b-', x, np.zeros(len(x)), 'k-', np.zeros(len(x)), 10*x, 'k-')


plt.annotate('', xy=(4.0,0), xytext=(3.9, 0),
            arrowprops=dict(facecolor='black', shrink=0.5),
            )
plt.text(3.4, 1, r'x')

plt.xlim(-pi/2, 0)
plt.ylim(-2,2)
plt.grid(True)

plt.show()
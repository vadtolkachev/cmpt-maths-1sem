import numpy as np

x_list = np.linspace(0, 1, 20000)
x_list = x_list[1:]

partial_amount = 0
for i in range(len(x_list)):
    partial_amount += 1/(1 + x_list[i]**2)

Integral = partial_amount * 0.0001

print(f'Integral = {np.round(Integral, 4)}')
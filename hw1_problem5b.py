import numpy as np
import matplotlib.pyplot as plt

delta = np.logspace(-16, 0, 17)

for x in [np.pi, 10**6]:

    original = np.cos(x + delta) - np.cos(x)

    transformed = -2 * np.sin((2*x + delta)/2) * np.sin(delta/2)

    difference = transformed - original

    plt.semilogx(delta, difference, 'o-', label=f'x = {x:g}')

plt.axhline(0, color='black', linewidth=0.8)

plt.xlabel(r'$\delta$')
plt.ylabel('Difference')
plt.title('Difference between the two')
plt.legend()
plt.grid(True)

plt.show()
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('2000.txt')

x = data[:,0]
y = data[:,1]

plt.plot(x,y)
plt.xlabel('$\lambda$')
plt.ylabel('Coeficiente de transmisión')

plt.show()

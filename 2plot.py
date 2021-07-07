import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

data1 = np.loadtxt('500.txt')
data2 = np.loadtxt('1000.txt')
data3 = np.loadtxt('2000.txt')
grill = []

grill.append(data1[:,1])
grill.append(data2[:,1])
grill.append(data3[:,1])


yl = [500,1000,2000]
xl = [0.1,0.3,0.5,1,5,10]

sb.heatmap(grill,xticklabels=xl,yticklabels=yl)

plt.show()

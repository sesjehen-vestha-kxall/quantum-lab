from model import hydrogen as hy

import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rd

print(hy.probability(hy.hydrogen_1s, r_b=0.68E-10))

SAMPLES = 200

rho, theta, phi = hy.linear_interval_values(samples=SAMPLES)

heatmap = np.zeros((SAMPLES,SAMPLES), dtype='float32')
for i in range(SAMPLES):
    for j in range(SAMPLES):
        heatmap[i,j] = hy.hydrogen_1s(rho[i])

plt.imshow(heatmap, cmap='hot', interpolation='nearest')
plt.show()
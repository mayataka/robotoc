import numpy as np
import matplotlib.pyplot as plt


def smoothl1(x, step_length, step_height, eps):
    assert eps > 0
    alpha = np.sqrt((0.5*step_length)**2 + eps)
    beta = step_height / (alpha - np.sqrt(eps))
    return beta * (alpha - np.sqrt((x-0.5*step_length)**2 + eps))


# step_length = 1.0
step_length = 0.12
step_height = 0.08
eps = 0.1
# eps = 1.0

x = np.arange(0, step_length, 0.001)
plt.plot(x, smoothl1(x, step_length, step_height, eps))
print(smoothl1(x, step_length, step_height, eps))
plt.show()
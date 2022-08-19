import numpy as np

X = np.array([[1, 2, 3],
              [2, 3, 4],
              [3, 4, 5],
              [4, 5, 6]])

w = np.array([10, 20, 30])

Y = X 
for i in range(X.shape[0]):
    Y[i, :] = X[i, :] + w
    
print(Y[1,:])

z = [1, 2, 3, 4, 5]
print(z[:])
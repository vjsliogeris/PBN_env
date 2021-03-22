"""utils.py
"""
import numpy as np


def integerize(X):
    """Turn a list of booleans in to an integer representation.
    """
    n_in = X.shape[0]
    s = 0
    for i in range(n_in):
        s += int(X[n_in-i-1]) * 2**i
    return s

def booleanize(X, l):
    """Turn an int to a boolean string of length l
    """
    output = np.zeros((l), dtype=bool)
    for i in range(l):
        h = 2**(l-i-1)
        if X >= h:
            X -= h
            output[i] = 1
    return output


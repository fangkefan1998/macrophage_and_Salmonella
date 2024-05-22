def normalize(array):
    import numpy as np

    max = array.max()
    min = array.min()
    normal_out = np.zeros([array.shape[0], array.shape[1]])

    for l in range(0, array.shape[0]):
        for w in range(0, array.shape[1]):
            normal_out[l,w] = (array[l,w] - min)/(max + 3000 - min)

    return normal_out
import os
import h5py


CHUNK_SIZE = 1000
MATRIX_DATASET = "X"


class Matrix(object):

    def __init__(self, file_name, dtype):
        self.file_name = os.path.abspath(file_name)
        self.dtype = dtype

    def add(self, X):
        X = X.astype(self.dtype)
        with h5py.File(self.file_name, "a") as hf:
            self.chunks = (CHUNK_SIZE, X.shape[1])
            if MATRIX_DATASET not in hf.keys():
                hf.create_dataset(MATRIX_DATASET, data=X, maxshape=(None, X.shape[1]), chunks=self.chunks)
            else:
                hf[MATRIX_DATASET].resize((hf[MATRIX_DATASET].shape[0] + X.shape[0]), axis=0)
                hf[MATRIX_DATASET][-X.shape[0]:] = X

    def shape(self):
        with h5py.File(self.file_name, "r") as hf:
            return hf[MATRIX_DATASET].shape

    def read(self):
        n = self.shape()[0]
        with h5py.File(self.file_name, "r") as hf:
            n = hf[MATRIX_DATASET].shape[1]
            for i in range(0, n, CHUNK_SIZE):
                chunk = slice(i, i + chunk_size)
                yield X[MATRIX_DATASET][chunk]

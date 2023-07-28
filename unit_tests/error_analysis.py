from typing import Any
import numpy as np
import matplotlib.pyplot as plt
import h5py


class ErrorAnalysis:
    def __init__(self, key):
        self.key = key

        path = "output/"
        self.org_data = h5py.File(f"{path}orgResult_{key}.h5", "r")["data"]
        self.new_data = h5py.File(f"{path}newResult_{key}.h5", "r")["data"]

        self.shape = self.org_data.shape
        assert np.all(self.new_data.shape == self.shape)


    def __call__(self):
        
        diff = np.abs( self.new_data[:] - self.org_data[:] )

        max_err = np.nanmax(diff)
        min_err = np.nanmin(diff)
        avg_err = np.nanmean(diff)

        self.diagnostics = [avg_err, min_err, avg_err]

        return avg_err
    
    def epicrisis(self):

        print("\n" + self.key)
        print("-"*len(self.key))
        print(f"> average error: {self.diagnostics[0]:8.4e}")
        print(f"> minimal error: {self.diagnostics[1]:8.4e}")
        print(f"> maximal error: {self.diagnostics[2]:8.4e}")

        print("")
    



##################### Loop Corrections #####################


loopCorr = ErrorAnalysis("LoopCorrection")
loopCorr()
loopCorr.epicrisis()


##################### Particle-Mesh Projections #####################

proj = ErrorAnalysis("PartMeshProjection")
proj()
proj.epicrisis()


##################### Fast Fourier Transforms #####################

fft = ErrorAnalysis("FasterFourierTransform")
fft()
fft.epicrisis()
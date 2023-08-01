from typing import Any
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import h5py

class ErrorAnalysis:
    def __init__(self, key):
        self.key = key

        path = "output/"
        try:
            self.org_data = h5py.File(f"{path}orgResult_{key}.h5", "r")["data"]
        except OSError:
            print("Referance data was not found. Run program without OpenMP first!")
            print("Exiting ...")
            exit()
        try:
            self.new_data = h5py.File(f"{path}newResult_{key}.h5", "r")["data"]
        except OSError:
            print("Test data was not found. Run program with OpenMP MPI!")
            print("Exiting ...")
            exit()

        self.shape = self.org_data.shape
        assert np.all(self.new_data.shape == self.shape)


    def __call__(self):
        tol = 1e-5
        msg = f"{self.key}: Values are dangerously close to zero"
        test = np.mean(np.abs(self.org_data)) > tol
        if not test:
            print(msg)
        
        diff = np.abs( self.new_data[:] - self.org_data[:] )

        max_err = np.nanmax(diff)
        min_err = np.nanmin(diff)
        avg_err = np.nanmean(diff)

        self.diagnostics = [avg_err, min_err, max_err]

        return avg_err
    
    def epicrisis(self):

        print("\n" + self.key)
        print("-"*len(self.key))
        print(f"> average error: {self.diagnostics[0]:8.4e}")
        print(f"> minimal error: {self.diagnostics[1]:8.4e}")
        print(f"> maximal error: {self.diagnostics[2]:8.4e}")

        print("")

    def plot(self):

        fig, ax = plt.subplots()
        data = np.mean(self.org_data, axis=2)
        # data = np.where(np.abs(data)>1e-3, data, np.nan)
        im = ax.imshow(data, cmap=cm.cool)#, norm=colors.LogNorm(vmin=1e-3))
        fig.colorbar(im)
        plt.show()
    



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


##################### Particle Velocity and Position Updates #####################

part = ErrorAnalysis("ParticleUpdate")
part()
part.epicrisis()

part.plot()
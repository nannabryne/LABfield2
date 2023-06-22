# LABfield2
Boilerplate examples and unit tests for investigating performance of LATfield2.

---

# TODO:
- [ ] Find a simple way of saving the parameters (`dim`, `npts` etc.) that was used to get the output?
- [ ] Find appropriate simulation parameters and parallel grid for automatically running tests

---

# Unit tests
We provide a number of unit tests for evaluating the performance (enhancement) of different parts of `LATfield2`. Each test represent a boilerplate example of some computation, placed in folders in `benchmarks/`. The results of said computation are to be saved (`fresh_output.h5`) and compared with the corresponding results (`org_output.h5`). The latter contains results produced with `LATfield2` _before_ any changes are made.[^1]

[^1]: I.e. prior to 26/06/23.

>Efficiency for a computation with $n$ compute processes and execution time $T$ is $$ \mathcal{E}(n)=\frac{n_{\mathrm{ref}} T_{\mathrm{ref}}}{nT},$$ where $n_{\mathrm{ref}}$ is the number of processes used in the original simulation and $T_{\mathrm{ref}}$ is the original simulation time. 


## **(a)** Field manipulation
> look-up folder: `benchmarks/field_manipulation/`


## **(b)** Execution of fast Fourier transform
> look-up folder: `benchmarks/fft_execution/`


## **(c)** Particle-mesh projection
> look-up folder: `benchmarks/particle_mesh_projection/`

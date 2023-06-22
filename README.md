# LABfield2
Boilerplate examples and unit tests for investigating performance of $\textup{LATfield2}$.

---

># TODO:
>- [x] Find a simple way of saving the parameters (`dim`, `npts` etc.) that was used to get the output?
>- [x] Find appropriate simulation parameters and parallel grid for automatically running tests (these can quite easily be changed)
>- [ ] Focus on the timing of each separate process (function call), rather than the runtime for the whole example
>- [ ] Deal with the issue that the reference data (`benchmarks/**/org_output.h5`) are way too large for git 
>   - maybe simply save samples in a tarball? 
>   - be ok that they will be saved locally?
>   - only save & compare slices of (or coarse-gridded) cubes?


---

# Unit tests
We provide a number of unit tests for evaluating the performance (enhancement) of different parts of $\textup{LATfield2}$. Each test represent a boilerplate example of some computation, placed in folders in `benchmarks/`. The results of said computation are to be saved (`fresh_output.h5`) and compared with the corresponding results (`org_output.h5`). The latter contains results produced with $\textup{LATfield2}$ _before_ any changes are made.[^1]

[^1]: I.e. prior to 26/06/23.

Note that each example (unit) includes several function calls (computations) to $\textup{LATfield2}$. A part of the result is saved and compared to the original result to make sure the changes do not alter the calculations.

>Efficiency for a computation with $n$ compute processes and execution time $T$ is $$ \mathcal{E}(n)=\frac{n_{\mathrm{ref}} T_{\mathrm{ref}}}{nT},$$ where $n_{\mathrm{ref}}$ is the number of processes used in the original simulation and $T_{\mathrm{ref}}$ is the original simulation time. 


## **(a)** Field manipulation
> look-up folder: `benchmarks/field_manipulation/`

The boilerplate example is that of a linear combination between two fields. 


## **(b)** Execution of fast Fourier transform
> look-up folder: `benchmarks/fft_execution/`



## **(c)** Particle-mesh projection
> look-up folder: `benchmarks/particle_mesh_projection/`

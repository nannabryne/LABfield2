# $\textup{LABfield2}$

Example codes to control the results from running code in [$\textup{LATfield2}$](https://github.com/nannabryne/LATfield2/tree/master) that is OpenMP parallelised.

(list of functions that are changed?)

## How to run
Suppose you are in the `unit_tests` directory.
To run the scripts using OpenMP, type
```bash
make; make run
```
and output will be written to `.h5` files in `output/`. 

If changes are made to the code in terms of number of particles, lattice resolution, etc., one may need to create a new result file to compare with. To run the scripts without OpenMP, one can simply type
```bash
make o=1; make run
```
and the original output files will be overwritten.

To run diagnostics, see `error_analysis.py`. One can also simply type
```bash
make analysis
```
and terminal print will be provided.

### Details
Flags are given be given to the compiler to specify what tests should be run.
- **`-DLOOPCORR`** $\leadsto$ `loopCorrectionSimple(...)` (see header `test_LoopCorrection.hpp`)
    >In `LATfield2/LATfield2_Lattice[_decl].hpp`, we test
    >1.  `Lattice::for_each(...)`
    > ---

- **`-DPROJECTION`** $\leadsto$ `partMeshProjectionSimple(...)` (see header `test_PartMeshProjection.hpp`)
    > In `LATfield2/particles/projections.hpp`, we test
    > 1. `scalarProjectionCIC_project(...)`
    > 2. `vectorProjectionCICNGP_project(...)`
    > 3. `symtensorProjectionCICNGP_project(...)`
    > ---
- **`-DFOURIER`** $\leadsto$ `fasterFourierTransformSimple(...)` (see header `test_FasterFourierTransform.hpp`)
    > In `LATfield2/LATfield2_PlanFFT_decl.hpp`, we test
    > 1. `PlanFFT::execute(FFT_FORWARD)`
    > 2. `PlanFFT::execute(FFT_BACKWARD)`
    > ---
- **`-DPART`** $\leadsto$ `particleUpdateSimple(...)` (see header `test_ParticleUpdate.hpp`)
    > In `LATfield2/particles/LATfield2_Particles.hpp`, we test
    > 1. `Particles::updateVel(...)`
    > 2. `Particles::moveParticles(...)`
    > ---


<!-- 

---
---
---
---
---
# $\textup{LABfield2}$
Boilerplate examples and unit tests for investigating performance of $\textup{LATfield2}$.

---

># TODO:
>- [x] Find a simple way of saving the parameters (`dim`, `npts` etc.) that was used to get the output?
>- [x] Find appropriate simulation parameters and parallel grid for automatically running tests (these can quite easily be changed)
>- [ ] Focus on the timing of each separate process (function call), rather than the runtime for the whole example
>- [ ] Deal with the issue that the reference data (`unit_tests/**/org_output.h5`) are way too large for git 
>   - maybe simply save samples in a tarball? 
>   - be ok that they will be saved locally?
>   - only save & compare slices of (or coarse-gridded) cubes? (I think this is a good option)
>- [ ] Mind the type of computer that is used -- find a way to be consistent about this
>- [ ] Hybrid programming: Masteronly approach (vector mode) ? or tasking (funneled or multithreaded)?


---

# Unit tests
We provide a number of unit tests for evaluating the performance (enhancement) of different parts of $\textup{LATfield2}$. Each test represent a boilerplate example of some computation, placed in folders in `unit_tests/`. The results of said computation are to be saved (`fresh_output.h5`) and compared with the corresponding results (`org_output.h5`). The latter contains results produced with $\textup{LATfield2}$ _before_ any changes are made.[^1]

[^1]: I.e. prior to 26/06/23.

Note that each example (unit) includes several function calls (computations) to $\textup{LATfield2}$. A part of the result is saved and compared to the original result to make sure the changes do not alter the calculations.

>Efficiency for a computation with $n$ compute processes and execution time $T$ is $\mathcal{E}(n)=\frac{n_{\mathrm{ref}} T_{\mathrm{ref}}}{nT}$, where $n_{\mathrm{ref}}$ is the number of processes used in the original simulation and $T_{\mathrm{ref}}$ is the original simulation time. We primarily use $n_{\mathrm{ref}}=64$. 


## **(a)** Field manipulation
> look-up folder: `unit_tests/field_manipulation/`

The boilerplate example is that of a linear combination between two fields. 


## **(b)** Execution of fast Fourier transform
> look-up folder: `unit_tests/fft_execution/`



## **(c)** Particle-mesh projection
> look-up folder: `unit_tests/particle_mesh_projection/`


# notes:
- Remember to think of race conditions!! Should maybe our code check whether we are in the danger zone or not?
- Performance enhancement & reduction in mem.  footprint -->

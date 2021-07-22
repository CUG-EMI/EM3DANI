# EM3DANI

EM3DANI is a package for isotropic/anisotropic 3D forward modeling of frequency-domain electromagnetic (CSEM and MT) data written in the [Julia language](http://julialang.org).

*  Authors: [Ronghua Peng](https://github.com/prhjiajie) and [Bo Han](https://github.com/hanbo1735) (China University of Geosciences (Wuhan)).

## License

The main part of EM3DANI is licensed under the [Apache Licence 2.0](http://www.apache.org/licenses/LICENSE-2.0), while the third-party [Dipole1D](https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/index.html) code which has been integrated into EM3DANI is licensed under the [GNU General Public License](http://www.gnu.org/licenses/).


## File structure
* **./doc** :        an instruction of file format.

* **./examples** :   contains subdirectories corresponding to different types of synthetic examples, including all those presented in the manuscript.

* **./src** :        source code.
* **./test** :       contains scripts for unit tests.

## Installation of Julia
EM3DANI is compatible with Julia v0.7 and later versions. We recommend to install v1.0.5, the long-term support (LTS) release.
### Windows systems
Go to the [Julia download page](https://julialang.org/downloads/) to download the Windows command line version (.exe) and install it.

### Linux systems
Although Julia is a cross-platform language, we strongly recommend to run EM3DANI under Linux rather than Windows. This is because some of the third-party packages utilized by EM3DANI such as **Dipole1D** and **MUMPS.jl** are more straightforward to complie under Linux.

There are three ways to install Julia on Linux:

* **Using precompiled binaries (recommended)**.	Go to the [Julia download page](https://julialang.org/downloads/) to download the generic Linux binaries (.tar.gz file).
	Then make sure that the Julia executable is visible for your system. To do this, first extract the .tar.gz file to a folder on your computer.
	Then you can either add Julia’s bin folder to your system PATH environment variable, or create a symbolic link to julia inside a folder which
	is on your system PATH, for example, by using the following command:

  `sudo ln -s <where you extracted the julia archive>/bin/julia /usr/local/bin/julia`


* **Compiling from source**. Assume that Git has been installed already, then we can grab the Julia sources from GitHub by using the following command:

  `git clone git://github.com/JuliaLang/julia.git`

  This will download the Julia source code into a julia directory in the current folder. The Julia building process needs the GNU compilation tools g++, gfortran, and m4, so make sure that you have installed them. Now go to the Julia folder and start the compilation process as follows:

  `cd julia`

  `make`


* **Using PPA for Ubuntu Linux**. Particularly, for Ubuntu systems (Version 12.04 or later), there is a Personal Package Archive (PPA) for Julia
	that makes the installation painless. All you need to do to get the stable version is to issue the following commands in a terminal session:

  `sudo add-apt-repository ppa:staticfloat/juliareleases`

  `sudo add-apt-repository ppa:staticfloat/julia-deps`

  `sudo apt-get update`

  `sudo apt-get install julia`

After a successful installation, Julia can be started by double-clicking the Julia executable (on Windows) or typing `julia` from the command line (on Linux). Following is the Julia's command line environment (the so-called REPL):


```jl
   _       _ _(_)_     |
  (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
  _  _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.0.5 (2019-09-09)
 _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
|__/                   |  

julia>
```

## Running the EM3DANI code
### Setting up the package environment
EM3DANI depends on several external packages (the so-called dependencies) which are not shipped with the package .zip file. These dependencies can be automatically resolved by activating and instantiating the package environment through [Julia's package manager (Pkg)](https://julialang.github.io/Pkg.jl/v1/). Go to the package directory, for example:  
`cd home/username/code/EM3DANI`

, and then enter the Julia REPL. Then press `]` from the Julia REPL you will enter the Pkg REPL which looks like
```jl
(v1.0) pkg>
```

, indicating that you are currently in the environment named v1.0, Julia 1.0's default environment. To switch to the package environment, just `activate` the current directory:
```jl
(v1.0) pkg> activate .
```
you will get:
```jl
(EM3DANI) pkg>
```
indicating that you are in the environment EM3DANI. The environment will not be well-configured until you `instantiate` it:
```jl
(EM3DANI) pkg> instantiate
```
. By doing so the dependencies listed in `Project.toml` and `Manifest.toml` can be automatically downloaded and installed. However, you may need to build one of those dependencies manually (see below for details).

The main direct dependencies of EM3DANI are three linear solver packages, namely **KrylovMethods.jl**, **MUMPS.jl** and **Pardiso.jl**. Both **KrylovMethods.jl** and **Pardiso.jl** are registered Julia packages. **MUMPS.jl** was registered, but it's not any more and has been renamed to [MUMPSjInv.jl](https://github.com/JuliaInv/MUMPSjInv.jl), and the reference of **MUMPS.jl** in `Manifest.toml` points to [our forked repository of the formerly registered version](https://github.com/CUG-EMI/KrylovMethods.jl).

 If you wish to use MUMPS as the linear solver, then you need to build **MUMPS.jl** manually after instantiating the package environment. First, find out where the Julia packages locate. By default, on Linux they are at (for example) `/home/username/.julia/packages/`. Then go to the MUMPS source folder (for example) `/home/username/.julia/packages/MUMPS/xxxxx/src`, you can find that there are two or three complier options files named like `compiler_options.in` or `compiler_options_XXX.in`. There are two options to `Make`: 
 
 * If you have **Intel compiler** (icc, ifort) combined with the **MKL library**, then simply type `Make` from the shell command line;
 
 * Otherwise, you need to have **GNU compiler** (gcc, gfortran) combined with the [OpenBLAS](http://www.openblas.net/) library preinstalled. Besides, **cmake** may be required. It seems that GNU compiler and cmake are usually missing in a freshly installed Linux system such as Ubuntu. They can be installed by typing the following commands from shell:
  
   `sudo apt-get install gcc`
   
   `sudo apt-get install gfortran`
   
   `sudo apt-get install cmake`
   
   The installation of OpenBLAS can be found [here](https://github.com/xianyi/OpenBLAS/wiki/Installation-Guide). Please track the installation directory, because the exact location of the complied OpenBLAS library must be given in `compiler_options_OpenBLAS.in`. For example,
   
   `LIBBLAS = /opt/OpenBLAS/lib/libopenblas_haswellp-r0.3.10.a`
   
   points to the complied library of OpenBLAS v0.3.10. After adjusting `compiler_options_OpenBLAS.in`, rename it to `compiler_options.in` and type `Make` from the shell command line.

To get back to Julia REPL from Pkg REPL, press `backspace` or `^C`.


### Testing and running the code
* First, you need to build the **Dipole1D** library. Go to the subdirectory of the **EM1DUtils** module, for example:  
`cd home/username/code/EM3DANI/src/EM1DUtils/deps`

  and then enter the Julia REPL, and "include" the script *build.jl*, which is like:
  ```jl
  julia> include("build.jl")
  ```

* Second, you need to let the EM3DANI package to be "loaded" by the current Julia environment. This is done by adding the parent directory of the package directory to  `LOAD_PATH`, a global environment variable of Julia. For example, the EM3DANI package is placed at `home/username/codes` on Linux or at `D:\\code` on Windows, then type the following command from the Julia REPL:

  ```jl
  julia> push!(LOAD_PATH,"/home/username/code")
  ```

  on Linux, or

  ```jl
  julia> push!(LOAD_PATH,"D:\\code")
  ```

  on Windows.   

* Third, go to the directory `.../EM3DANI/test`, you will see several Julia scripts which are for unit tests, run the one named *runtest.jl* by typing the following command from the Julia REPL:
  ```jl
  julia> include("runtest.jl")
  ```
  you will get the following information if no error occurs:
  ```jl
  Test Summary: | Pass Total
  EM3DANI       |  149   149
  ```

* Finally, go to the directory where the running script loated, and run the script by typing the following command (for example) from the Julia REPL:

  ```jl
  julia> include("runFwd.jl")
  ```

### Writing a running script
For each numerical example contained in the directory `./examples`, one or more running scripts named `runFwd.jl` or `runFwd_parallel.jl` have been provided. These scripts are well documented. A user can modify them to get his/her own.
Following are explanations of several parameter settings which are not that easy to understand.

#### Reference model
The forward modeling algorithm is based on a secondary field formulation, thus a (1D) background/reference model must be provided to perform forward modeling. This is done by defining a variable of the composite type **RefModel**. It has two fields: *sig1D* and *depth1D*, containing the conductivities (S/m) and top depths (m) of an isotropic layered model, respectively. For example, the following codes define a variable of type **RefModel** named *refModel*, representing a marine 1D model with three layers (air-water-sediment):

```
refModel = RefModel(zeros(3), zeros(3))    # construct a variable of type RefModel
refModel.sig1D   = [1e-8, 3.3, 1.0]        # the conductivities of air, sea water, and sediment
refModel.depth1D = [-100000, 0, 1000.0]    # the top depths of air, sea water, and sediment layers
```
Because the top boundary of the reference model must coincide with that of the total model domain (read from the model file), the top depth of the first layer can be set to `-emMesh.origin[3]` for correctness and convenience:
```
refModel.depth1D = [-emMesh.origin[3], 0, 1000.0]
```
where **emMesh** is a variable of composite type **EMTensorMesh**, containing model & mesh parameters reading from the model file.


#### Linear solver
The parameters of linear solver are stored in a variable of composite type **DirectSolverParm** or **IterativeSolverParm**. For example, the following codes define a parameter set for the MUMPS direct solver:

```
lsParm = DirectSolverParm()
lsParm.solverName = :mumps    # can be :mumps or :mklpardiso.
lsParm.sym        = 1         # 0=unsymmetric, 1=symm. pos def, 2=general symmetric
lsParm.ooc        = 0         # 0=in-core, 1=out-of-core
lsParm.saveFac    = false     # whether or not saving the decomposation factors
```
for solving a symmetric, positive definite linear system with "In-Core" mode, and not saving the matrix decomposation factors. Because the above parameter set is exactly the default value of a **DirectSolverParm** *constructor*, the above codes are equivalent to:
```
lsParm = DirectSolverParm()
```
But if you want to use MKL Pardiso instead of MUMPS, then you have to set the solver name explicitly:
```
lsParm = DirectSolverParm()
lsParm.solverName = :mklpardiso
```
It is worth mentioning that the memory requirement of a direct solver can be dramatically reduced
by enabling the out-of-core (OOC) mode (setting lsParm.ooc to 1), in which the matrix factors are stored in disk rather than RAM of the computer. 

With regarding to the number of threads used by a direct solver (shared memory parallelism), although the type DirectSolverParm has a field *nThread* to specify it for MKL Pardiso, it is recommend to specify it by setting the global environment variable *OMP_NUM_THREADS* or *MKL_NUM_THREADS*, for example:
```
ENV["OMP_NUM_THREADS"] = 4
ENV["MKL_NUM_THREADS"] = 4
```


For iterative solvers, on the other hand, the following is an example
```
lsParm = IterativeSolverParm()
lsParm.iterMethod = :qmr        # can be :bicgstb or :qmr.
lsParm.tol        = 1e-7        # relative residual tolerance
lsParm.prec       = :aphi       # preconditioning method, only :aphi is allowed
lsParm.maxIter    = 1000        # maximal iteration number
```
which defines a QMR solver preconditioned with *A-phi*, with a relative residual tolerance of 1e-7 and a maximal iteration number of 1000.

### Parallel execution of the code
Two different levels of parallelism have been implemented in EM3DANI: 

* **Solving the large linear system using direct solvers:** at a low level the direct solver MUMPS or MKL Pardiso can solve the linear system of equations in parallel by exploiting a multi-threaded linear algebra library such as OpenBLAS and MKL BLAS (*shared memory parallelism*). 
  
  To utilize this type of parallelizaion, besides choosing a direct solver as the linear solver, you need to do nothing else. The code will involve multiple threads in forward computation automaticlly. In our experience, the direct solver will involve all available threads in solving equations if there is no explicit restriction on the number of threads, which can not always obtain the highest efficiency due to the relatively poor scalability of direct solvers (see [Puzyrev et al., 2016](https://www.sciencedirect.com/science/article/pii/S0098300416300164) [Han et al., 2018](https://library.seg.org/doi/10.1190/geo2017-0515.1) and among others) as well as the possible interferences among threads. The optimal number of threads generally increases with the problem size and the capacity of the computing platform. In our presented numerical examples, for example, the optimal number of threads is between 4 and 8 for the MT 1D anisotropic example (with a grid size of 22\*40\*102) and between 8 and 16 for the CSEM 1D example (with a grid size of 52\*114\*67).
 
  The number of threads is stored in the global environment variable *OMP_NUM_THREADS* or *MKL_NUM_THREADS* (when the parallelization is offerred by the MKL BLAS library), which can be set from the shell like
  
   `shell> export OMP_NUM_THREADS = 4`

   or from the Julia REPL like

   `julia> ENV["OMP_NUM_THREADS"] = 4`

* **Computing for multiple frequencies & transmitters:** at a higher level forward computations of multiple frequencies & transmitters can be distributed to different computer processors (*distributed parallelism*) by using Julia’s parallel computing mechanism.

  To utilize this type of parallelizaion, first you need to launch multiple worker processes by either starting Julia like

  `shell> julia -p 8`

  or adding processes within Julia (recommended) like

  `julia> addprocs(8)`
  
  Then you need to call the parallel forward modeling function **parsolveEM3DFwd** instead of the sequential one **solveEM3DFwd** (please refer to the various scripts named `runFwd_parallel.jl` within the `examples` directory).

  There are two schemes for distributing the multiple forward computations: in terms of only frequencies (regardless of tranmitters/polarization modes) (Scheme 1), and in terms of frequencies and tranmitters (Scheme 2). The code is able to choose the optimal one automatically: if a) a direct solver is used as the linear solver or b) the number of tranmitters/polarization modes is less than two (which is impossible for MT) or c) the number of frequencies is divisible by the number of worker processes, Scheme 1 will be chosen, otherwise Scheme 2 will be chosen.

In principle, the two levels of parallelism are independent of each other. However, they can affect each other in practice due the limited computing resources. For example, if we employ *n* worker processes to do the forward computation, each one solves the linear system using a direct solver with *m* threads involved, then there may be *n\*m* threads working simultaneously and the instantaneous peak memory usage may reach *n* times that of a single worker process, which can resulting unexpected poor scablity if the computer is not powerful enough.

Therefore, to obtain a fine-tune parallel performance, a user need to balance the number of worker processes and the number of threads used by the direct solver (if chosen), taking the problem size, number of frequencies and tranmitters/polarization modes, and computer capacity into consideration.

### Using the Julia plotting tools
There are a few Julia plotting scripts in each numerical example directory (subdirectory of `.../EM3DANI/examples`) for plotting the data. These are developed by using the Julia package [Plots.jl](http://docs.juliaplots.org/latest/) and one of its "backends" [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl), both of which have been set as dependencies of EM3DANI (see the project file `Project.toml`). Thus, once EM3DANI is "instantiated", these plotting scripts can be directly run to generate figures without having to install any extra package. Nevertheless, We suggest one install [Juno](https://junolab.org/)/[Atom](https://atom.io/), so that a figure can be explicitly displayed in an interactive panel.



## Notes
* EM3DANI is not supposed to deal with complex topography since it employs a rectangular mesh.

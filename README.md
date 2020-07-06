# EM3DANI

EM3DANI is a package for isotropic/anisotropic 3D forward modeling of frequency-domain electromagnetic (CSEM and MT) data written in the [Julia language](http://julialang.org).

*  Authors: Ronghua Peng and Bo Han (China University of Geosciences (Wuhan)).

## License

The main part of EM3DANI is licensed under the [Apache Licence 2.0](http://www.apache.org/licenses/LICENSE-2.0), while the third-party [Dipole1D](https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/index.html) code which has been integrated into EM3DANI is licensed under the [GNU General Public License](http://www.gnu.org/licenses/).


## File structure
* **./doc** :        an instruction of file format.

* **./examples** :   contains subdirectories corresponding to different types of synthetic examples, including all those presented in the manuscript. Each subdirectory has its own README.

* **./src** :        source code.

## Setting up the Julia environment
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

After a successful installation, Julia can be started by double-clicking the Julia executable (on Windows) or typing the `julia` from the command line (on Linux). Following is the Julia's command line environment (the so-called REPL):


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

## Installing prerequisite Julia packages

EM3DANI utilizes three third-party Julia packages as the linear solver, namely **KrylovMethods.jl**, **MUMPS.jl** and **Pardiso.jl**. Although you may wish to use only one of them, all the three packages need to be "added" to your Julia environment, which is done through Julia's package manager (Pkg). To Enter the Pkg REPL, press `]` from the Julia REPL and you should see a similar prompt:
```jl
(v1.0) pkg>
```
To get back to the Julia REPL, press backspace or ^C.

* **KrylovMethods.jl** is a *registered* Julia package, thus it can be added simply by typing `add KrylovMethods` from the **Pkg REPL**. However, the current registered version does not contain an implementation of **QMR** method. Therefore, we recommend to add our forked version by typing:

 `add https://github.com/CUG-EMI/KrylovMethods.jl`

* **MUMPS.jl** was a *registered* Julia package, but it's not any more and has been renamed as [MUMPSjInv.jl](https://github.com/JuliaInv/MUMPSjInv.jl). Therefore, we recommend to add our forked version by typing:

 `add https://github.com/CUG-EMI/MUMPS.jl`

 If you wish to use MUMPS as the linear solver, then you need to build **MUMPS.jl** manually after adding it. First, find out where the Julia packages locate. By default, on Linux they are at (for example) `/home/username/.julia/packages/`. Then go to the MUMPS source folder (for example) `/home/username/.julia/packages/MUMPS/xxxxx/src`, you can find that there are two or three complier options files named like `compiler_options.in` or `compiler_options_XXX.in`. There are two options to `Make`: If you have *Intel compiler (icc, ifort) combined with the MKL library*, then simply type `Make` from the shell command line; otherwise, you need to have *GNU compiler (gcc, gfortran) combined with the [OpenBLAS](http://www.openblas.net/) library* preinstalled (the [installation of OpenBLAS](https://github.com/xianyi/OpenBLAS/wiki/Installation-Guide) is quite straightforward), and rename the corresponding complier options file `compiler_options_OpenBLAS.in` as `compiler_options.in` before typing `Make`.

* **Pardiso.jl** is a *registered* Julia package, thus it can be added simply by typing `add Pardiso` from the Pkg REPL, and it will be built automatically. Currently we have only utilized the MKL version of Pardiso, thus you need to have the MKL library preinstalled if you want to use
Pardiso as the linear solver while running EM3DANI (but it does not affect the building process).


## Running the EM3DANI code

* First, you need to let the EM3DANI package to be "loaded" by the current Julia environment. This is done by adding the parent directory of the package directory to  `LOAD_PATH`, a global environment variable of Julia. For example, the EM3DANI package is placed at `home/username/codes` on Linux or at `D:\\code` on Windows, then type the following command from the Julia REPL:

 `push!(LOAD_PATH,"/home/username/code")`   (on Linux)

 or

 `push!(LOAD_PATH,"D:\\code") `        (on Windows)   


* Second, go to the directory where the running script loated, and run the script by typing the following command (for example) from the Julia REPL:

 `include("runFwd.jl")`

There are several example running scripts in subdirectories of the directory `./example`, which are well documented. Please
refer to them for how to call the functions of EM3DANI to perform forward modeling.

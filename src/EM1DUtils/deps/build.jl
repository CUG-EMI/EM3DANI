#------------------------------------------------------------------------------
#
# script to build dipole1d shared library
# (c) PRH, CUG-IGG, 24 July, 2015
#
#------------------------------------------------------------------------------
using BinDeps

@BinDeps.setup

# setup filepath
depsDir = splitdir(Base.source_path())[1]
srcDir  = joinpath(depsDir, "src")
libDir  = joinpath(depsDir, "lib")

printstyled("=== Building Dipole1D shared library === \n", color=:blue)
println("depsDir = $(depsDir)")
println("srcDir  = $(srcDir)")
println("libDir  = $(libDir)")

if !isdir(libDir)
    printstyled("=== Creating library directory === \n", color=:blue)
    mkdir(libDir)
    if !isdir(libDir)
        error("Creating library directory failed.")
    end
end

if Sys.isunix() # is_unix()
    src01 = joinpath(srcDir, "FilterModules.f90")
    src02 = joinpath(srcDir, "Dipole1D.f90")
    src03 = joinpath(srcDir, "callDipole1D.f90")
    dirfile = joinpath(libDir, "Dipole1D")

    #@build_steps begin
        printstyled("=== Fortran compilier Info ===\n", color=:blue)
        run(`gfortran --version`)
        run(`gfortran -O2 -march=native -fPIC -shared $(src01) $(src02) $(src03) -o $(dirfile)`)

    #end

end

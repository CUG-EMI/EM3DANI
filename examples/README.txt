There are overall five numerical examples, including an MT 1D isotropic example, an MT 3D isotropic example, an MT 1D anisotropic example, a CSEM 1D anisotropic example, and a CSEM 3D anisotropic example. The latter three are those presented in the manuscript.

Files:
*.dat       data file
*.mod       model file
*.resp      response file
*.jl        Julia script to run forward modeling


****************** 1 - MT 1D isotropic example ******************
A three-layer isotropic model, the result can be easily validated by an analytic solution.


****************** 2 - MT 3D isotropic example ******************
The COMMEMI-3D2 model. Various solutions of forward modeling of this model can be found in the literature, for example:

  1) Mackie, R. L., T. R. Madden, and P. E. Wannamaker, 1993, Three-dimensional magnetotelluric modeling using difference equations — Theory
  and comparisons to integral equation solutions: Geophysics, 58, 215–226, doi: 10.1190/1.1443407.
  
  2) Siripunvaraporn, W., G. Egbert, and Y. Lenbury, 2002, Numerical accuracy of magnetotelluric modeling: A comparison of finite
  difference approximations: Earth Planets Space, 54, 721–725, doi: 10.1186/BF03351724.
  
Thus, our result can be verified.

The folder 'plotFwdResp' contains Julia scripts to read and plot the output response.


****************** 3 - MT 1D anisotropic example ******************
Corresponds to the model shown in Figure 3 in the manuscript. The two Julia scripts 'runFwd.jl' and 'runFwd_parallel.jl', are used to run forward modeling in sequential and parallel modes, respectively.

The folder 'AnalyticalSolution' contains analytical solutions, and the folder 'plotCompareNumAna' contains Julia plotting scripts to reproduce Figure 4 of the manuscript.


****************** 4 - CSEM 1D anisotropic example ******************
Corresponds to the model shown in Figure 5 in the manuscript. It contains two cases: VTI and dipping anisotropy. The script 'runFwd_parallel_mumps.jl' is used to run forward modeling in parallel mode, using MUMPS as the solver, while 'runFwd_parallel_qmr.jl' is also used to run forward modeling in parallel mode, but using QMR as the solver.

The folder './VTI/plotCompareNumAna' contains Julia plotting scripts to reproduce Figure 6 a, b, c and d of the manuscript.

The folder './Dip_30/plotCompareMFV_vs_FE' contains Julia plotting scripts to reproduce Figure 6 e and f of the manuscript.


****************** 5 - CSEM 3D anisotropic example ******************
Corresponds to the model shown in Figure 7 in the manuscript. It contains two cases: Azimuthal and dipping anisotropy. Each case contains multiple specific models but only one running script. Therefore, in order to obtain the forward response of every model, you must adjust the running script repeatedly (modify the name of model file and output file).

The folder './Dipping/compareDipVary' contains Julia plotting scripts to reproduce Figure 8 a, b, d and e of the manuscript.

The folder './Azimuthal/compareStrikeVary' contains Julia plotting scripts to reproduce Figure 10 a, b, d and e of the manuscript.

The folders './Dipping/plotFieldVector' and './Azimuthal/plotFieldVector' are for reproducing Figure 9 and 11. To do this, the grid electric fields at grid sampling points must be output as a file (*.field) (See the last a few lines in the running script). (Because this type of file generally has very large size, we do not provide our results).

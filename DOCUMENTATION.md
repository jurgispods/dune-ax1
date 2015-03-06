# Documentation of configuration options for simulations
## Configuration at run-time or compile-time
Most of the configurations and settings for one of the axon simulations are made in a config file.
The reason for this is that the compilation of an executable takes a considerable amount of time (on the order of a few minutes), such that configuration at run-rime is desirable. 

However, some of the settings (e.g., the number of ion species) have to be made at compile-time, in the files [`dune/ax1/common/constants.hh`](dune/ax1/common/constants.hh) or in one of the configuration headers residing in `dune/ax1/acme2_cyl/configuration`. This README gives a short summary about the most important settings that can be made in either a config file (run-time parameters) or in the config headers (compile-time parameters).

## Config headers (compile-time)
### constants.hh
This section goes through the file constants.hh, which contains global constants and settings for all of the executables and different configurations.

#### BaseGrid
```
#if USE_GRID==1
typedef Dune::YaspGrid<2,Dune::TensorProductCoordinates<double,2> > BaseGrid;
#elif USE_GRID==2
typedef Dune::UGGrid<2> BaseGrid;
#endif
```
Depending on the preprocessor define `USE_GRID` (defined in `src/Makefile.am`), the type of the "base grid" is chosen. All recent code relies on a `YaspGrid`, the `UGGrid` used in older versions of the code is not supported anymore.

#### NUMBER_OF_SPECIES
```
//! \brief Number of ion species
static const unsigned int NUMBER_OF_SPECIES = 3;
```
This global constant gives the number of concentrations to be used. In my productive simulations, I have always used three (Na+, K+, Cl-) species. For validations or convergence tests, switching to a single species has proven useful as well.


#### DO_MULTIDOMAIN_BLOCKING
```
#ifdef MULTIPLE_MEMBRANE_ELEMENTS
#define DO_MULTIDOMAIN_BLOCKING (0)
#else
#define DO_MULTIDOMAIN_BLOCKING (1)
#endif
```
The default (and highly recommended!) default value is true. When using only one membrane element layer, the number of DOFs at each grid vertex is the same (i.e. `NUMBER_OF_SPECIES + 1`, one unknown per each concentration variable plus one potential unknown). This fact allows blocking together all DOFs at a given vertex into block matrices of size `NUMBER_OF_SPECIES + 1` x `NUMBER_OF_SPECIES + 1`. Such a blocking is highly efficient in conjunction with preconditioners like ILUn, which allow to exactly invert those small block matrices.

The preprocessor define `MULTIPLE_MEMBRANE_ELEMENTS` switches off multidomain blocking. This is necessary since with more than one membrane element layer, the inner membrane vertices each have only one potential DOF associated with them, but not concentration unknowns. A limitation in PDELab currentlt only allows for blocking when all sub-block matrices have equal size, which is not the case with multiple membrane elememts. Consequently, blocking has to be switched off, and solving very large problems is not possible, since the solver will not converge in reasonable time.

#### AX1_BLOCKSIZE
```
//! \brief Blocksize in ISTL Backends
#if DO_MULTIDOMAIN_BLOCKING==1
// Use this for a fully blocked matrix
static const unsigned int AX1_BLOCKSIZE = NUMBER_OF_SPECIES+1;
#else
// Use this for a non-blocked matrix
static const unsigned int AX1_BLOCKSIZE = 1;
#endif
```
The blocksize is directly depending on the blocking flag and should not be changed manually.

#### USE_CYLINDER_COORDINATES
```
//! \brief Flag for using cylinder coordinates in acme2_cyl
#ifdef USE_CYLINDER_COORDS
static const bool USE_CYLINDER_COORDINATES = true;
#else
static const bool USE_CYLINDER_COORDINATES = false;
#endif
```
One nice thing about the current implementation is that one can change between a plain 2D setup or a cylinder-symmetrical setup by virtue of a preprocessor define `USE_CYLINDER_COORDS`, which sets or unsets the corresponding global bool flag. We make use of this in the `src/Makefile.am` in order to generate two executables from the same source code, one for the plain 2D setup (`acme2_2d`) and one for the cylinder-symmetrical setup (`acme2_cyl`).

#### ION_NAMES
```
//! \brief enum for ion species
enum { Na = 0, K = 1, Cl = 2, Ca = 3 };
static const std::vector<std::string> ION_NAMES = { "Na", "K", "Cl"};
```
These helper structures allow for meaningful indexing and printing of the ion species.

#### SUBDOMAIN_NAMES
```
//! \brief enum for element subdomains (formerly known as 'element groups')
enum { CYTOSOL = 0, ES = 1, MEMBRANE = 2};
static const std::vector<std::string> SUBDOMAIN_NAMES = { "CYTOSOL", "ES", "MEMBRANE"};
```
Analoguus to the ion species, these constructs allow for meaningful indexing and printing of the element subdomains. Note that these element groups to not correspond to the actual grid subdomains, as cytosol and extracellular space are combined in an electrolyte grid subdomain (see below).

#### MD_SUBDOMAINS
```
//! \brief enum for subdomains (when using dune-multidomaingrid)
/* Important: Domains 0,..,nSubdomains-1 must correspond to the subdomain indices returned
 * by subGV.grid().domain(); -1 is a special value for the complete domain; further custom
 * values can be added with decreasing values starting with -2. This is especially useful
 * for identifying special gridfunctions which live on exotic subsets of entities (like my
 * beloved membrane interface GFs)
 */
enum GridDomains { DOMAIN_ALL = -1, DOMAIN_ELEC = 0, DOMAIN_MEMB = 1, DOMAIN_MEMB_INTERFACE = -2};
static const std::vector<std::string> MD_SUBDOMAINS = { "ELEC", "MEMB"};
```
As a complication to the previous element groups, there also exist grid subdomains. The only two subdomains that actually exist in the `MultidomainGrid` are electrolyte and membrane. But it turned out to be useful to create virtual subdomains (like "all" or "membrane interface") for the output of gridfunctions.


#### Physical constants
```
// physical and mathematical constants in SI units
//! pi === 3
static const double con_pi   = 3.1415926535897932384626433;
//! elementary charge e in [C = A s]
static const double con_e    = 1.602176487e-19;
//! Boltzmann constant k_B = R/N_A in [J / K]
static const double con_k    = 1.380650e-23;
//! vacuum permittivity eps0 in [C / ( V m )]
static const double con_eps0 = 8.854187817e-12;
//! Avogadro constant N_A in [1/mol]
static const double con_mol  = 6.02214129e23;
```
These common physical constants are used frequently throughout the code.

#### debug streams
```
//! \brief DUNE debug streams for convenient debug outputs
static const Dune::DebugLevel APPL_MINLEVEL = 1;
Ax1ParallelDebugStream<Dune::DebugStream<1, APPL_MINLEVEL> > debug_verb(std::cout);
Ax1ParallelDebugStream<Dune::DebugStream<2, APPL_MINLEVEL> > debug_info(std::cout);
Ax1ParallelDebugStream<Dune::DebugStream<3, APPL_MINLEVEL> > debug_warn(std::cerr);
Ax1ParallelDebugStream<Dune::DebugStream<3, APPL_MINLEVEL> > debug_minimal(std::cout);

enum LogLevel { VERB=0, INFO=1, WARN=2, MINIMAL=3};

// Set level to 2 to deactivate Jochen's output
static const Dune::DebugLevel JOCHEN_MINLEVEL = 1;
Ax1ParallelDebugStream<Dune::DebugStream<1, JOCHEN_MINLEVEL> > debug_jochen(std::cout);
```
We created several different output streams, which can be separately activated or deactivated (set first template parameter below `*MINLEVEL`). Alternatively, one can increase the `*MINLEVEL` to only let streams with first template parameter >= `*MINLEVEL` generate output.


### configuration headers
The configuration headers reside in `dune/ax1/acme2_cyl/configurations`, each configuration in a separate subfolder. The available configurations are:

- default: PNP equations including all subdomains (cytosol, membrane, ES)
- ES: PNP equations on a single extracelullar domain
- laplace: VC (volume conductor) equation on a grid with or without membrane
- mori: EN (electroneutral) equations as developed by Yoichiro Mori (see e.g. http://dx.doi.org/10.2140/camcos.2009.4.85) on a grid with or without membrane
- myelin: Deprecated; originally created for simulation of myelinated axons, but this functionality has been fully integrated into the default config
- step, test_scales: Test and validations configurations with only a single ion species on a single grid domain

Each of the configurations has at least three important header files: one `*_config.hh` configuration header, and parameter headers for concentrations (`*_nernst_planck_params.hh`) and potential (`*_poisson_boundary.hh`).

#### *_config.hh
This file contains some compile-time flags and typedefs specifig to the given configuration. Some of those flags are just default values that can be overridden by run-time flags from a config file. Notable are:

##### numberOfSpecies
The number of ion species this configurations works with. Must match with `NUMBER_OF_SPECIES` from `constants.hh`, otherwise it will result in an error.

##### con_diffWaterNa, con_diffWaterK, con_diffWaterCl
Diffusion coefficients for each of the ion species

##### Traits
Traits struct containing some important typedef like the type of the concentrations parameter class (`NERNST_PLANCK_PARAMETERS`) or the local operator for the electrolyte subdomain (`ELEC_OPERATOR`). Note that these typedefs might not be fully specified according to the template parameters, but leave some of its template parameters open. This is possible by means of the new "alias declaration" functionality from C++11 (see http://stackoverflow.com/a/2795024/2429404, for example) and the `using` keyword. This comes handy when considering the complex type instantiation chain in multiple setup files and avoids forward declarations.

Two important flags from the `Traits` struct are:

- useMembraneContributions: Whether or not the membrane elememts contribute to the finite element matrix/residual. It has proven useful to deactivate membrane contributions in setup like VC ("laplace") or EN ("mori), where a membrane is included in the grid for geometrical reasons, but the equations are only defined on the electrolyte subdomains
- useImplicitMembraneFlux: Switch between explicit and implicit membrane flux handling. This is a compile-time flag for now, since both versions are implemented in different classes with different types, so a recompilation is necessary when switching from implicit to explicit or vice versa.

#### *_nernst_planck_params.hh / *_poisson_boundary.hh
These two parameter classes contain all equation-specific coefficients of the PDEs (with respect to a very general convection-diffusion equation) as well as boundary condition types and values. They are the core classes for the different models considered (PNP, VC, EN).

## Config files (run-time)
A number of different config files can be found in the `src/config*/` subfolders. They contain an overwhelming number of different parameters, most of which are reasonable default values which do not have to be modified. In the following, I will try to give an overview over all available parameters and their purpose.

### [general] section
The `[general]` contains most of the simulation parameters. It would greatly profit from a restructuring into subsections, but I never found time to do this.

#### Input/Output
##### Basic Output
```
outputDir = out-acme2_cyl_par_mori_test
```
This is the (relative or absolute) output directory to which all of the simulation output will be written to. It must be created if it does not exist, otherwise the simulation will exit with an error.

```
vtkOutput = no
vtkSubSamplingLevel = 0
```
These parameters control whether VTK (www.vtk.org) files are written out. The subsampling level lets you output more than the actual unknowns by virtual refinement of the grid.

```
gnuplotOutput = yes
fullGnuplotOutput = no
```
Gnuplot output plots out plain-text files that can be visualized by Gnuplot. Very handy for debug purposes, as they can be view "on-the-fly", while a simulation is running. Full gnuplot output should be avoided for large grids or in production runs, as the amount of generated data is considerably larger.

```
hdf5Output = yes
fullHdf5Output = yes
```
HDF5 (www.hdfgroup.org/HDF5) is a space-efficient binary file format (some call it a file system within a single file) and the output method of choice for large parallel runs. The amount of output can be reduced by switching off full HDF5 output.

```
oneTimeGridfunctionOutput = yes
```
One-time gridfunction output is done once at time=0. It includes all time-independent data like grid geometry, parallel grid paritioning, or membrane groups (in the case of a myelinated axon, this includes the location of nodes of Ranvier and myelin sheath as well as the membrane permittivity).

```
diagnosticsFile = diagnostics.dat
```
The diagnostic file contains various diagnostic output like time, time step size, number of Newton iterations, defect of the finite element residual before and after one time step, etc..

```
outputTimeInterval = 0.0
```
This value specifies the time interval between two subsequent outputs. By setting it to zero, you require output at every time step.

##### Loading/saving simulation states
```
# This flag simply loads a state
loadState = no
```
If you want to load a state file containing the solution from a previous simulation, this flag needs to be set to true.

```
doExtrapolateXData = no
```
A special case is when you want to load from a grid that is not matching the current grid. The current implementation allows to load from a grid with the same points in y-direction, but with different resolution in x-direction. This is handy when loading an equilibrium state, since in this case the equilibration simulation can be carried out on a much coarser grid (one cell in x-direction is sufficient, as all the dynamics happen in y-direction, perpendicular to the membrane orientation). In this case, set this flag to true and the coarse solution will be extrapolated onto the fine grid.

```
loadGridFileName = simulation_states/equilibrium/cyl/yasp_square10mm.dgf
```
This file specifies the grid from which to load in form of a DGF (Dune Grid Format) file.

```
continueSimulation = no
```
This flag tells to continue the history of the last simulation, i.e. existing data files in the output dir are not cleared, new data is simply appended; if false previous history is "forgotten".

```
loadFilename = simulation_states/equilibrium/cyl/yasp_square10mm.dat
```
This file contains the actual state file containing the solution to load.

```
loadFromParallelRun = no
parallelRunNumberProcessors = 64

```
These parameters are only needed when loading a parallel simulation state, consisting of `parallelRunNumberProcessors` state files, into a sequential simulation. Of course, the macro grids of the previous parallel run and the current sequential run must match.

```
saveState = yes
saveFilename = simulation_states/hackepeter/acme2_cyl_par_mori_test.dat
```
These parameters control whether or not and where the final simulation result should be saved.

```
doCheckpointing = yes
checkpointInterval = 10
```
Checkpointing is useful for long-running simulations. If activated, it will produce complete snapshots of the current solution at the specified interval, which can be used in case of a crash to load the intermediate solution back into a new simulation.

```
createGridFile = yes
saveGridFileName = acme2_cyl_par_mori_test_equilibrium.dgf
```
These flags specify if a DGF grid file should be written out, which is required when loading the obtained solution onto a new grid.

```
loadChannelStates = no
```
When loading from a previous simulation, this flag controls whether not only the concentration/potential DOFs, but also the internal channel states ("gating particles") should be loaded from the state file. 
Please note: In the EN model by Y. Mori, additional state variables exist on the membrane; for these, the load/store functionality has not yet been implemented!

##### Parallel output
```
rootOutputNode = 0
outputRootNodeOnly = yes
```
When running in parallel, we usually choose rank 0 to be the output node, but this can be changed. Furthermore, it is possible to activate output on all ranks by setting the `outputRootNodeOnly` flag to false.

##### Numerical debug output
```
printMatrix = no
printRhs = no
printResidual = no
```
These parameters allow to write out the finite element matrix and right-hand side in each Newton iteration. Additionally, the residual at every grid vertex can be printed out after every time step. All of these outputs are quite space- and time-consuming and should be switched off in a productive run.

#### General numerical and grid parameters
```
useRowNormPreconditioner = no
```
The row norm preconditioner normalized the finite element matrix and rhs by dividing each row by the maximum entry. This is a rudimentary, but sometimes necessary preconditioning step in addition to the actual precinditioner/solver combination. 
One situation where it is necessary is when using different boundary condition types for concentrations and potential at the same location.

```
configName = default
```
This speciefies the static (compile-time) configuration to use. Can be one of {default, ES, laplace, mori, myelin, step, test_scales}.

##### Grid generation
The parameters associated with grid generation are not easy to describe, as they are not independent of each other. For a complete and reliable description of the grid generation process, please take a look at the class `Ax1GridGenerator`.

```
refineXDirection = no
refineYDirection = yes
```
These flags work together with the first command line argument `level`, which speciefies the level of refinement of the coarse grid. Refinement in x- and y-direction can be switched on or off separately. However, depending on the other parameters `refineYDirectionGeometrically` and `useMembrane` described below, the above parameters may work differently.

```
refineYDirectionGeometrically = yes
```
This flag is very important when using cylinder coordinates, as we want a fine resolution close to the membrane, but a much coarser grid towards the outer boundary. The geometric refinemnt does just that.

```
doAdditionalGlobalRefine = no
```
Setting this to true will result in a global refinement after the coarse grid has been created. In contrast to the usual `grid.globalRefine(int)` mechanism from the Dune grid interface, this one excludes the membrane in y-direction. Therefore, the number of membrane element layers do not change, but all the other elements are refined as usual.

```
# radial cylinder refinement compensates for increasing volumes in y-direction
useRadialCylinderRefinement = no
# Start coordinate of radial refinement
radialCylinderRefinementStart = 15.0e3
# Shall there be an interval of smooth geometric refinement between membrane and
# the beginning of the radial cylinder refinement (radialCylinderRefinementStart)?
transitionIntervalGeometricRefinement = yes
# If not, use nElementsTransitionInterval equidistant elements instead
nElementsTransitionInterval = 20;
# Desired max. ratio between two extracellular cell volumes
max_volume_ratio = 1.e3
# Max. number of points to use in radialCylinderRefinement
max_radial_points = 100
```
The "radial refinement" mechanism was a poor workaround to deal with the increasing grid element volumes in radial direction when using cylinder coordinates. The preferred way to deal with this is by using volume scaling (see below), so radial refinement can be seen as obsolete.
The flag `useRadialCylinderRefinement` allows to keep the volumes of elements constant (or within certain bounds) in y-direction. Of course, this come at the price that the mesh sizes get smaller and smaller, so it is not possible use this mechanism with a large value of `ymax`. 
The coordinate at which to start radial refinement (`radialCylinderRefinementStart`) can be specified as well.
Between the membrane and this coordinate, one may choose a smooth refinement strategy (flag `transitionIntervalGeometricRefinement`) or, alternatively, the number od equidistant elements in between (`nElementsTransitionInterval`). 
The parameter `max_volume_ratio` controls the maximum allowed ratio between two element volumes over the radial refinement range. But a harder condition is given by `max_radial_points`, which is the maximum allowed number of points spent on the radial refinement range, even if this violates the `max_volume_ratio` condition.

##### Membrane parameters
```
useMembrane = yes
```
This flag specifies whether we consider a geometry with a membrane or not.

```
numberMembranes = 1
```
The number of membranes in this setup. Their location is specified by `yMemb`, see below.

```
refineMembrane = no
nMembraneElements = 1
``` 
Whether or not to refine the membrane(s) and if yes, how many membrane element layers do we want? Note: The default case (no refinement, 1 membrane element layer) has the great advantage of enabling matrix blocking (see above), and one membrane element layer has proven to be enough for accuracy purposes, so we strongly advice to use the default values.

##### Various numerical parameters
```
mV_output = yes
```
This flag enables conversion of dimensionless potential to units of mV in output files.

```
gammaPot = 1.0
gammaCon = 1.0
```
These are deprecated damping parameters used in a previously implemented operator-splitting approach. It is not used in the fully-coupled Newton approach, which is the preferred numerical scheme.

```
# one of {step, squarePulse, trianglePulse, smoothStep}
#initConc = default
```
Deprecated parameter to choose from different initialization procedures. By now, the initialization class used is specified in the configuration header's `Traits` struct at compile-time.

##### Time step control
```
adaptiveTimeStep = yes
#timeStepFile = simulation_states/load_mori/timesteps.dat
#timeStepEvery = 10
timeStepLowerLimitNewtonIt = 10
timeStepUpperLimitNewtonIt = 30
```
The flag `adaptiveTimeStep` controls whether or not to use adaptive time-stepping. 
If set to false, one may specify a time step file from which to read in the predefined time steps to use. An additional parameter `timeStepEvery` tells the program to only use every `timeStepEvery`-th time step value from the specified file. If not file is specified, a constant time step size is used throughout the simulation, as handed over by the second command line argument.
When using an adaptive time step, the time step size in increased or decreased depending on the number of newton iterations needed in the last time step. The last two parameters control the thresholds for increasing (#iterations < `timeStepLowerLimitNewtonIt`) or decreasing (#iterations >= `timeStepUpperLimitNewtonIt`) the time step size.

##### Solver thresholds and numerical tuning parameters
```
newtonReduction = 1e-5
newtonAbsLimit = 1e-2
newtonMaxIt = 10
```
Three parameters are involved in the Newton convergence criterion: Either the current defect of the finitel element residual must be smaller than `newtonAbsLimit`, or the achieved reduction (ratio between current and initial defect of finite element residual) smaller than `newtonReduction`. In any case, the maximum allowed number of iterations is given by the hard threshold `newtonMaxIt`. If the convergence criterion was not satisfied withon the allowed number of iterations, and the maximum number of restarts (see parameter `maxNumberNewtonRestarts`) has been exhausted, the simulation exits with an error.

```
newtonMinLinReduction = 1e-5
```
The Newton implementation from PDELab usually picks a reasonable convergence criterion for the linear solver used internally. However, it has proven useful to require a minimum reduction from the linear solver, since a premature convergence might lead to non-convergence in the following Newton iteration.

```
maxNumberNewtonRestarts = 0
```
For some problems, the Newton iteration might not converge because the time step was chosen too large. By allowing a certain number of Newton restarts, we can retry the non-converged iteration using half of the previous time step size.

```
linearSolverMaxIt = 500
```
The maximum number of linear solver iterations in each Newton iteration.

```
disableLineSearch = no
```
This flag allows to disable the line search strategy in the Newton iteration

```
doReorderMatrix = no
```
This flag should be switched on when using multiple membrane element layers. As carried out earlier, using more than one membrane element layer prevents matrix blocking, but the matrix entries can still be ordered by grid vertices way, which results in a favorable band-diagonal structure.

```
doVolumeScaling = yes
#volumeScalingYThreshold = -1
```
Volume scaling assigns a scaling factor to each DOF, compensating for increasing cell volumes algebraically. This is used as a replacement for radial refinement. It allows to use large radial domain size without losing numerical stability. By default, a reasonable extracellular starting coordinate is chosen, which is about one axon radius from the membrane. It can be adjusted through the parameter `volumeScalingYThreshold`.

```
#potentialResidualScalingFactor = 1e6
```
This parameter scales the potential residual by a certain factor. It can be used to put emphasis on the importance of the potential variables compared to the concentration variables, otherwise the default value (1.0) is recommended.

#### Domain parameters
```
cellWidth = 500.
y_memb = 500.
d_memb = 5.
xmax = 10.e6
ymax = 10.e6
ymin = 0.
```
These parameters control the extent of the computational domain in x-direction (`xmax`) and y-direction(`ymax`) as well as the axonal cytosol radius (`cellWidth`) and the membrane thickness (`d_memb`). 
By means of `ymin`, the y-component of the grid origin can chosen different from zero, which is interesting for purely extracellular simulations.
One important parameter is `y_memb`, which is actually a vector of length `nMembranes`, containing the lower staring coordinates of each membrane in the domain. Of course, its values should be consistent with `ymin` and `cellWidth`.

#### Mesh sizes
```
dx = 100.e3
# dy only takes effect when logarithmic grid is chosen, otherwise dx is determined by the level
dy = 100.e3
dy_min = 0.5
n_dy_min = 10
dy_cell = 100.
dy_cell_min = 100.
n_dy_cell_min = 0
```
By default, the mesh sizes are equidistant in x-direction (`dx`). When using the default grid creation procedure with geometric refinement in y-direction, `dy` and `dy_min` give the maximum and minimum extracellular mesh sizes and `dy_cell` and `dy_cell_min` the maximum and minimum intracellular mesh sizes in y-direction. 
The parameters `n_dy_min`and `n_dy_cell_min` specify how many equidistant elements of the minimum sizes will be placed directly adjacent to the extra- or intracellular membrane interfaces, respectively.

#### Various simulation setup parameters
```
closedCell = yes
```
This depracated flag enforced Neumann-0 boundary conditions on all intracellular boundaries. Since the `[boundary]` section now allows to explicitly specify the boundary condition typed separately for each boundary, this flag is obsolete.

```
# Initialize channels to the current membrane potential
automaticChannelInitialization = yes
# Value to use if no automatic initialization is desired
channelRestingPotential = -65.
```
This flag will initialize all the ion channels to their steady-state values with respective to the initial membrane potential. If it is deactivated, a fixed membrane potential of value `channelRestingPotential` will be used instead.

##### EN model flags
One important note about the EN model: There exist two different executables for this model, `acme2_cyl_par` together with the `mori` configuration uses the fully-coupled Newton approach, `acme2_cyl_par_mori` uses the operator-splitting scheme ((in the latter case, the only allowed configuration is `mori`).

```
useMoriChargeLayerContribution = yes
```
This flag activated the EN model. Must be true when using the `acme2_cyl_par_mori` executable. Must be true when using the `acme2_cyl_par` executable together with the `mori` configuration 

```
useMoriOperatorSplit = no
```
Chooses operator-splitting as the numerical scheme; if deactivated, it is solved fully-coupled by Newton iteration. Must be true when using the `acme2_cyl_par_mori` executable. Must be false when using the `acme2_cyl_par` executable.

```
disableMoriFluxCalculation = no
```
When using the EN model, this flag deactivates the additional ion flux contributions stemming from the additional membrane state variables in the EN model. Useful when loading previously generated fluxes as boundary conditions; in this case we can safely deactivate the calculations of fluxes, as they are not used, to speed up the computation.

```
useElectroneutralityTerminationCriterion = no
```
Mori suggests to use an electroneutrality condition as the termination criterion for the operator split. This can be switched off; in this case, the parameter `newtonMaxIt` is abused to specify the number of iterations for the operator-split iteration. I have used this for most of my calculations, since I suspected a bug in my implementation of the electroneutrality check.

```
useMoriCorrection = no
```
Mori also suggests a postprocessing of the concentrations/potential in order to recover the exponential Debye layer profile. In a previous version, I applied this postprocessing after every solver step and used the obtained values for the solution of the next iteration. This is, however, a serious intrusion and might make the system ill-posed or prevent convergence. Therefore, this flag is depracated in the current implementation. It should be changed to only be used for output, probably making the result comparable to the PNP solution.

```
removeDebyeLayerVertices = yes
minimalMembraneDistance = 709
fillRemovedRange = yes
```
Since the EN model does not require to resolve the Debye layer, I implemented a mechanism to first create a grid like in the PNP case, and afterwards removing the Debye layer vertices. This is useful in order to compare PNP and EN at the same points. 
If the flag `removeDebyeLayerVertices` is true, all the vertices from the extracellular membrane interface up to a distance of `minimalMembraneDistance` will be removed. 
If the flag `fillRemovedRange` is true, the removed range will be filled with a small number of relatively coarse elements, such that the mesh size does not vary too much due to the removed vertices.

```
doConcentrationPostprocessing = no
```
Mori suggests a postprocessing step to avoid global charge accumulation in his finite volume scheme. I did not observe any problems in my finite element code, but the postprocessing step can be switched on nevertheless.

```
# This factor scales the membrane flux. Only used for speedup of equilibration phase!
dummy_value = 1
```
This notoriously named parameter is a scaling factor that multiplies the membrane flux. It can be used to speed up equilibration simulations (and in effect it has only an effect if both `doEquilibration` is true and the current time is smaller than `tEquilibrium`, see below for both these parameters).

### [boundary] section
#### Membrane flux parameters
```
fullyImplicitMembraneFlux = yes
```
When true, membrane flux is calculated implicitly, explicitly otherwise. Note: The corresponding compile-time flag from the configuration headers must be consistent, otherwise the simulation will exit with an error!

```
verboseMembraneFluxClasses = no
```
When set to true, the membrane flux classes will print out debug information about the flux calculations for every membrane element. Should be deactivated in a productive run.

```
writeBoundaryOutput = yes
```
If set to true, the program will generate additional debug output about the boundary values, which is especially useful for debugging membrane flux calculations.

#### Boundary conditions
```
debyeLayerWidth = 0.
```
This parameter controls the width of the Debye layer in y-direction; important when using different boundary condition types on the Debye layer boundaries. In this case, the side boundaries consist of three parts (e.g. `leftCytosol`, `leftExtracellular`, `leftDebye` for the left boundary, see below). If set to zero (default), the side boundaries only consist of two parts (e.g. `leftCytosol` and `leftExtracellular` for the left boundary, see below).

```
useTimeDependentBoundaryValuesCon = no
useTimeDependentBoundaryValuesPot = no
```
These flags tell the program to load boundary values for concentrations and/or potential from a previously generated file (see below)

```
couplePotentialNeumannBoundaryToConcentration = yes
```
Tells program to calculate the potential Neumann boundary conditions as the sum of concentration boundary conditions; mandatory for EN model

```
loadBoundary = membrane
``` 
When loading either potential or concentration boundary conditions from a file, this parameter specified on which boundary the values should be loaded (can only be one boundary!)

```
membraneOneSided = no
```
Deprecated flag, imposes flux conditions only on the extracellular side of the membrane interface, should be false.

```
collapseToLineSource = no
originalRadius = 505.
```
This flag can be used in the VC setup to test the cylindrical axon vs an axon that was collapsed to a line source. In this case, the `originalRadius` parameter has to be provided with the location of the extracellular membrane interface. Since both setup give approx. the same results, it is recommended to keep the flag deactivated.

```
boundaryConLoadFilename = simulation_states/load_mori_complete/total_flux_na.dat simulation_states/load_mori_complete/total_flux_k.dat simulation_states/load_mori_complete/total_flux_cl.dat
#boundaryPotLoadFilename = simulation_states/load_laplace/total_flux_dx10e3.dat
```
These parameters specify the file names from which to load concentration or potential bondary values, respectively.

### [boundary.concentration] section
```
#membrane = Neumann
top = Dirichlet
bottom = Neumann
leftCytosol = Neumann
rightCytosol = Neumann
leftExtracellular = Neumann
rightExtracellular = Neumann
#leftDebye = Dirichlet
#rightDebye = Dirichlet
```
All these parameters control the boundary condition typed for the concentrations for each boundary location.

Important: Use row-norm preconditioner when using Dirichlet constraints for concentrations and Neumann constraints for potential on side boundaries, otherwise the (`NUMBER_SPECIES+1`x`NUMBER_SPECIES+1`) block matrices will be (close to) singular in ILU decomposition!


### [boundary.potential] section
```
membrane = Neumann
top = Dirichlet
bottom = Neumann
leftCytosol = Neumann
rightCytosol = Neumann
leftExtracellular = Neumann
rightExtracellular = Neumann
#leftDebye = Dirichlet
#rightDebye = Dirichlet
```
All these parameters control the boundary condition typed for the potential for each boundary location.

### [equilibration] section
```
doEquilibration = no
```
Whether or not to perform an equilibrium simulation (in this case, active ion channels are deactivated).

```
forceEquilibration = no
```
Normally, equilibration simulations are carried out with only one element in x-direction. This flag forces the equilibration calculation, even with more than one element in x-direction.

```
tEquilibrium = 1000.0
dtEquilibrium = 10.0
```
`tEquilibrium` speciefied the time up to which the equilibration is calculated. At this time point, the active channels are opened and a "normal" simulation is started. The used time step size up to this time is given by `dtEquilibrium`; tEquilibration must be a multiple of dtEquilibrium.
We recommend carrying out the equilibration simulation until tend = `tEquilibrium` and then saving the equilibrium state and then load this state in a new, non-equilibration simulation.

```
Na_leak = 0.13
K_leak = 0.87
```
Percentage of the available leak channel conductances from the total leak conductance `leak`; sum must be 1!

### [stimulation] section
```
stimulation = yes
```
Switches stimulation on or off.

```
memb_flux_stimulation = yes
position_x = 150.e3
position_y = 0.
```
When this flag is true, the stimulation happens by adding a certain amount to the membrane flux at the x-coordinate `position_x` of the first membrane; otherwise, sodium in injected directly into the domain at position (`position_x`,`position_y`).

```
t_inj_start = 10.0
t_inj_end = 2.0e3
```
Start and end times of the stimulation

```
I_inj = 1.1e13
```
Strength of the injection when using sodium injection (i.e., not using membrane flux injection).

```
memb_flux_ion = Cl
memb_flux_inj = 2e7
```
When using membrane flux injection, `memb_flux_ion` specifies the ion species for which to increase the flux, and `memb_flux_inj` gives the maximum value of flux stimulation, which is increased smoothly between `t_inj_start` and `t_inj_end`. A positive value means flow of the respective ion in positive y-direction (from CY to ES), i.e. a depolarizing current for Cl- and a hyperpolarizing current for Na+/K+.


### [solution_in] section
```
Na =  12.
K =  125.
Cl = 137. 
pot = -65.
```
Initial values for intracellular concentrations (units mM) and potential (units mV).

### [solution_ex] section
```
Na = 100.
K =  4.
Cl = 104.
pot = 0
```
Initial values for extracellular concentrations (units mM) and potential (units mV).

### [membrane] section
This `[membrane]` section allows to specify an arbitrary number of subsection (e.g. myelin, nodes_of_ranvier) in order to specify different channel densities on different parts of the axon. In this section, we specify some general or default parameters for the whole axon.

```
isActive = yes
```
Completely activate or deactivate ion channels; actually, this is a vector of size `nMembranes`!

```
useMultipleGroups = no
```
Whether or not to use multiple membrane groups

```
defaultGroup = myelin
```
Default membrane group; used for x intervals with unspecified membrane groups.


```
smoothPermittivities = yes
```
Smooth permittivity coefficient at membrane group interfaces; this only has an effect when membrane groups are resolved fine enough.

```
loadMultipleStates = yes
```
When using more than one membrane group, we can (and have to) load equilibrium states corresponding to each of the membrane groups.

```
interpolateMultipleStates = no
```
This flag indicates to interpolate the multiple loaded equilibrium states at membrane group transitions.

### [membrane.subgroup] sections
```
[membrane.node_of_ranvier]
start = 500.e3
width = 1.e3
stride = 1000.e3
```
Each subgroup defines `start`, `width` and `stride` parameters; in this example, the first node of Ranvier extends from 500e3 to 501e3, the second from 1500e3 to 1501e3, and so on.

```
loadFilename = simulation_states/hackepeter/acme2_cyl_par_y10e6_EQUILIBRIUM.dat
```
If we set the flag `loadMultipleStates` to true, we must provide a file name for each membrane subgroup.

```
permittivity = 2.
#d_memb = 500.
```
For each membrane subgroup, we can specify values for the `permittivity` and the membrane thickness `d_memb`. Note that the grid will not actually have varying membrane thicknesses, but it will rather calculate an effective permittivity, see my PhD thesis for details.

```
dx = 100.
```
Mesh size in x-direction for this subgroup.

```
smoothTransition = no
geometricRefinement = yes
dx_transition = 100.
n_transition = 10

```
`smoothTransition` says whether or not the mesh size will be varied smoothly towards the adjacent membrane group transitions. If yes, we can specify `geometricRefinement` towards the transition, which will use a  `n_transition` equidistant cells of mesh size `dx_transition`.

```
leak   = 0.5
Nav    = 120.
Kv     =  36.
```
Maximum conductances for the total `leak` conductance or various active channel conductances (units  are [mS/cm^2])



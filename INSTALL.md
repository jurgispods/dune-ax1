# Installation

## Dune dependencies

The minimum requirement for this Dune module are the Dune [core modules](http://www.dune-project.org/coremodules.html)

- dune-common
- dune-geometry
- dune-grid
- dune-istl
- dune-localfunctions

and the following additional modules:

- [dune-pdelab](http://www.dune-project.org/pdelab) and dune-typetree
- [dune-multidomaingrid](https://github.com/smuething/dune-multidomaingrid)
- [dune-multidomain](https://github.com/smuething/dune-multidomain)

The dune-ax1 module contains branches that correspond to the respective branches in dune-pdelab, i.e. the master branch works with the master branch of dune-pdelab, the branch `maint-pdelab_release_2.0` works with the 2.0 release branch of dune-pdelab (und corresponding release branches of the other modules), and so on. The latest tested state on the master branch is given by the following Git revisions:

- dune-common: commit 5b0bdc28589f4d66b8b56c6a4d781cf19f10774f
- dune-geometry: commit 3dafcc20186428e9669a90d69eae3bcc8b08d9a1
- dune-grid: commit 12091a82828433ab050c70c0a87e643197993f97
- dune-istl: commit 35b879c3f0d003dcbe743c7ca84f6142751c6b48
- dune-localfunctions: commit 25a1116ad55bc183ecd5d66232564c85c22cb455
- dune-typetree: commit 6730ea7fd6c44b6206555e664be071c40dc10b48
- dune-pdelab: commit 9f4d7a979c5369cfe2232d77e034253202f69a8b
- dune-multidomaingrid: commit 406c5c750650ebaf1b5811e9fb6012899b238064
- dune-multidomain: commit 2aa14eeb57d3f94563b926c637ebfe9235112c9d

Additionally, there are a set of patches in the stuff/ folder, which have to be applied on top of some of the Dune modules mentioned above.
Suppose you are at the root of your Dune tree, the following commands will apply the patch for dune-pdelab:

```
cd dune-pdelab
patch -p1 < ../dune-ax1/stuff/pdelab_ax1.patch
```

## Additional dependencies

- external libraries as required by the other Dune modules (may include gfortran, blas, boost)
- OpenMPI
- hdf5 (parallel, compatible to OpenMPI version)
- SuperLU (direct solver; not required, but useful)

## Compiling Dune stack from source

For instructions on how to install the Dune stack, see the INSTALL_DUNE.txt file or have a look at the [Dune website](http://www.dune-project.org/doc/installation-notes.html).
You will generally need a custom opts file, some examples can be found in the stuff/ folder.
Please note that this module uses the autotools build system, the new CMake system introduced to Dune recently is not supported.
The resulting call to `dunecontrol` for installation may look like:

```
dune-common/bin/dunecontrol --no-cmake --module=dune-ax1 --opts=my_opts_file.opts all
```

You will need a recent compiler like gcc 4.7, clang 3.2 or newer for this.


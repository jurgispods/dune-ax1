#!/bin/bash

#git clone http://git.dune-project.org/repositories/dune-common
#git clone http://git.dune-project.org/repositories/dune-geometry
#git clone http://git.dune-project.org/repositories/dune-grid
#git clone http://git.dune-project.org/repositories/dune-istl
#git clone http://git.dune-project.org/repositories/dune-localfunctions
#git clone http://git.dune-project.org/repositories/dune-pdelab
#git clone https://git.gitorious.org/dune-multidomain/dune-multidomain.git
#git clone https://git.gitorious.org/dune-multidomaingrid/dune-multidomaingrid.git

cd dune-common
git checkout 892d2043ab7554c69e6c99b985a33ced4d3bc627
cd ..

cd dune-geometry
git checkout 46256e328b40438dd27630464d9167aeb2b4a779
cd ..

cd dune-grid
git checkout 01fa5120dbaad9c9bfd76863b05603a44798bb25
cd ..

cd dune-istl
git checkout fcff302d4f9e0629047a8664fd72c793832f4e3e
cd ..

cd dune-localfunctions
git checkout 29416af5f482949f333fdb1c57653cca7e66821a
cd ..

cd dune-pdelab
git checkout f2e997429c5b103ef9b8684cb3a1ba26c0e01bf8
patch -p1 < ../dune-ax1/stuff/pdelab_ax1.patch
cd ..

cd dune-multidomain
git checkout a61c621fc9f9fa962e459ae93a5ae4eb42a47b20
cd ..

cd dune-multidomaingrid
git checkout 29a858b83fb65def59fa30d12d5a55e11606756e
cd ..

cd dune-ax1
git checkout next
cd ..

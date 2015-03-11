/*
 * acme2_cyl_output.hh
 *
 *  Created on: Aug 11, 2011
 *      Author: jpods
 */
#ifndef DUNE_AX1_ACME2CYL_OUTPUT_HH
#define DUNE_AX1_ACME2CYL_OUTPUT_HH


#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#if HAVE_HDF5
#include <dune/ax1/common/hdf5_tools.hh>
#endif
#include <dune/ax1/common/ax1_gridfunctions.hh>

#include <dune/ax1/acme2_cyl/common/acme2_cyl_simloader.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_solutionvectors.hh>
#include <dune/ax1/acme2_cyl/common/acme2_cyl_geometrytools.hh>



//template<class Acme2CylTraits, class GV, class SubGV, class GFS_CON, class GFS_POT, class U_CON, class U_POT, class PHYSICS>
template<class Acme2CylTraits>
class Acme2CylOutput
{
  public:

    // Boilerplate code to ease type usage
    typedef typename Acme2CylTraits::Coord Coord;
    typedef typename Acme2CylTraits::Real Real;
    typedef typename Acme2CylTraits::GridView GV;
    typedef typename Acme2CylTraits::SubGridView SubGV;
    typedef typename Acme2CylTraits::MultiGFS MultiGFS;
    typedef typename Acme2CylTraits::GFS_CON_SUB GFS_CON; // SubGFS!
    typedef typename Acme2CylTraits::GFS_POT_SUB GFS_POT; // SubGFS!
    typedef typename Acme2CylTraits::U U;
    //typedef typename Acme2CylTraits::U_CON U_CON;
    //typedef typename Acme2CylTraits::U_POT U_POT;
    typedef typename Acme2CylTraits::Physics PHYSICS;

    // Not used at the moment, but might be useful instead of using SelectComponentGridFunctionAdapter
    template<int SPECIES>
    struct GFS_ION
    {
      typedef Dune::PDELab::GridFunctionSubSpace<GFS_CON,Dune::PDELab::TypeTree::TreePath<SPECIES> > type;
    };

    const int SUBSAMPLING_POINTS;

    const std::string vtkOutputBase_Domain;
    const std::string vtkOutputBase_Membrane;

    // A bunch of typedef's for all the gridfunctions
    struct Acme2CylOutputTraits
    {
      typedef typename Acme2CylTraits::GFS_CON_SUB GFS_CON_SUB; // SubGFS!
      typedef typename Acme2CylTraits::GFS_POT_SUB GFS_POT_SUB; // SubGFS!

      // DGFs
      typedef Dune::PDELab::DiscreteGridFunction        <GFS_POT,U> DGF_POT;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U> DGF_POT_GRAD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,DGF_POT_GRAD,PHYSICS> DGF_POT_GRAD_ELEC;

      typedef Dune::PDELab::VectorDiscreteGridFunction  <GFS_CON,U> DGF_CON;
      typedef Dune::PDELab::Ax1VectorDiscreteGridFunctionGradient<GFS_CON,U> DGF_CON_GRAD;

      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON, PHYSICS> DGF_CON_MD;
      typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD, PHYSICS> DGF_CON_GRAD_MD;

      typedef ChargeDensityGridFunction<DGF_CON_MD,PHYSICS> GF_CD;
      typedef Ax1ConductivityGridFunction<typename Acme2CylTraits::PARAMETERS_POT,DGF_CON_MD,PHYSICS> GF_CONDUCTIVITY;

      // TODO Do not restrict and interpolate again, can this be simplified?
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_MD,1> DGF_SINGLE_CON_MD;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_GRAD_MD,GV::dimension> DGF_SINGLE_CON_GRAD_MD;

      // Meta grid function for the electrolyte ion flux taking only gridfunctions defined on the elec subdomain
      typedef IonFluxGridFunction<DGF_CON, DGF_CON_GRAD, DGF_POT_GRAD_ELEC, PHYSICS> GF_ELEC_FLUX;
      //typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ELEC_FLUX,2> GF_SINGLE_ELEC_FLUX;

      typedef typename Acme2CylTraits::GF_MEMB_FLUX BGF_MEMB_FLUX;
      typedef MembraneFluxGridFunction_DiffusionTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DIFF;
      typedef MembraneFluxGridFunction_DriftTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DRIFT;
      typedef MembraneFluxGridFunction_LeakFlux<BGF_MEMB_FLUX> BGF_MEMB_FLUX_LEAK;
      typedef MembraneFluxGridFunction_VoltageGatedFlux<BGF_MEMB_FLUX> BGF_MEMB_FLUX_VOLTAGE_GATED;

      typedef typename Acme2CylTraits::GF_MORI_FLUX BGF_MEMB_FLUX_MORI;
      typedef Ax1MembraneCurrentGridFunction<typename Acme2CylTraits::PARAMETERS_CON> BGF_MEMB_CURRENT;
      typedef Ax1ParametersToMembraneFluxAdapter<typename Acme2CylTraits::PARAMETERS_CON> BGF_PARAM_MEMB_FLUX;

      typedef Ax1ParametersToNeumannAdapter<typename Acme2CylTraits::PARAMETERS_POT,BGF_PARAM_MEMB_FLUX> BGF_PARAM_POT_NEUMANN;

      // Grid functions living on the membrane (only one component for the flux in normal direction)
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX,PHYSICS> GF_MEMB_FLUX;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DIFF,PHYSICS> GF_MEMB_FLUX_DIFF;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DRIFT,PHYSICS> GF_MEMB_FLUX_DRIFT;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_LEAK,PHYSICS> GF_MEMB_FLUX_LEAK;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_VOLTAGE_GATED,PHYSICS>
        GF_MEMB_FLUX_VOLTAGE_GATED;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_MORI,PHYSICS>
        GF_MEMB_FLUX_MORI;
      typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_CURRENT,PHYSICS> GF_MEMB_CURRENT;

      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX> GF_SINGLE_MEMB_FLUX;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX_DIFF> GF_SINGLE_MEMB_FLUX_DIFF;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX_DRIFT> GF_SINGLE_MEMB_FLUX_DRIFT;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX_LEAK> GF_SINGLE_MEMB_FLUX_LEAK;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX_VOLTAGE_GATED>
        GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_FLUX_MORI>
        GF_SINGLE_MEMB_FLUX_MORI;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_MEMB_CURRENT> GF_SINGLE_MEMB_CURRENT;
      typedef Dune::PDELab::IntersectionGridFunctionSelectComponentAdapter<BGF_PARAM_MEMB_FLUX> GF_SINGLE_PARAM_MEMB_FLUX;

      // Grid functions living on the whole domain (one component per coordinate axis)
      typedef Ax1ElectrolytePlusMembraneIntersectionToMultiDomainGridFunction<GV,PHYSICS,GF_ELEC_FLUX,BGF_MEMB_FLUX> GF_ION_FLUX;
      // Use SelectComponentGridFunctionAdapter, which is able to extract more than one component from the father GF!
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ION_FLUX,GV::dimension> GF_SINGLE_ION_FLUX;

      typedef MembranePotentialGridFunction<DGF_POT, PHYSICS> GF_MEMB_POT_MD;
      typedef Ax1MultiDomainGridFunctionRestriction<SubGV,GF_MEMB_POT_MD,PHYSICS> GF_MEMB_POT;
      typedef MembranePotentialGridFunction_MembGV<SubGV,DGF_POT,PHYSICS> GF_MEMB_POT_MEMBGV;

      typedef ChannelGridFunction<GV,Real,PHYSICS> GF_CHANNEL;

      typedef Ax1MembraneElementIndexGridFunction<GV,PHYSICS> GF_MEMB_GROUPS;

      // Analytical solutions
      typedef typename PHYSICS::Traits::SOLUTION_CON ANALYTICAL_SOLUTION_CON;
      typedef typename PHYSICS::Traits::SOLUTION_POT ANALYTICAL_SOLUTION_POT;
      typedef Dune::PDELab::SelectComponentGridFunctionAdapter<ANALYTICAL_SOLUTION_CON,1>
        SINGLE_ANALYTICAL_SOLUTION_CON;

      // Geometry helper gridfunctions
      typedef Acme2CylPlain2DAreaGridFunction<GV,Real> GF_PLAIN2D_AREA;
      typedef Acme2CylVolumeGridFunction<GV,Real> GF_VOLUME;
      typedef Acme2CylBaseSurfaceGridFunction<GV,Real> GF_BASE_SURFACE;
      typedef Acme2CylSideSurfaceGridFunction<GV,Real> GF_SIDE_SURFACE;
      typedef Acme2CylSurfaceGridFunction<GV,Real> GF_SURFACE;
      typedef Ax1PermittivityGridFunction<GV,PHYSICS> GF_PERMITTIVITY;
      typedef Ax1ProcessorGridFunction<GV> GF_PARTITION;

      // Boundary value functions
      typedef Ax1BoundaryValueFunction<DGF_POT,PHYSICS> GF_POT_BOUNDARY;
      typedef Ax1BoundaryValueFunction<DGF_CON_MD,PHYSICS> GF_CON_BOUNDARY;
    };


    typedef Acme2CylOutputTraits Traits;
    typedef typename Traits::DGF_POT::Traits::DomainType DomainType;

    struct DiagnosticInfo
    {
      DiagnosticInfo(bool hasAnalyticalSolution_)
      : time(0.0),
        iterations(0), dt(0.0),
        l2ErrorPot(0.0), maxErrorPot(0.0),
        l2ErrorCon(NUMBER_OF_SPECIES,0.0), maxErrorCon(NUMBER_OF_SPECIES,0.0),
        tEquilibrium(0.0),
        maxDiffCon(0),
        maxDiffPot(0),
        hasAnalyticalSolution(hasAnalyticalSolution_)
      {}

      enum ReductionFunction { Function_Identity = 0, Function_Min = 1, Function_Max = 2,
        Function_Mean = 3, Function_Sum = 4  };

      // diagnostic variables for one time step
      double time;
      int iterations;
      double dt;
      Real l2ErrorPot;
      Real maxErrorPot;
      std::vector<Real> l2ErrorCon;
      std::vector<Real> maxErrorCon;
      double tEquilibrium;

      std::vector<Real> maxDiffCon;
      std::vector<Real> maxDiffPot;
      bool hasAnalyticalSolution;

      std::map<std::string, Real> debugData;
      std::map<std::string, int> debugDataReductionFunctions;


      Real getl2ErrorCon(int i) const
      {
        return l2ErrorCon[i];
      }

      void setl2ErrorCon(int i, Real err)
      {
        l2ErrorCon[i] = err;
      }

      Real getMaxErrorCon(int i) const
      {
        return maxErrorCon[i];
      }

      void setMaxErrorCon(int i, Real err)
      {
        maxErrorCon[i] = err;
      }

      void setMaxDiffConPot(Real maxDiffCon_, Real maxDiffPot_)
      {
        maxDiffCon.push_back(maxDiffCon_);
        maxDiffPot.push_back(maxDiffPot_);
      }

      std::string getHeader()
      {
        std::stringstream header;
        int count = (hasAnalyticalSolution ? 6 : 4);
        if(hasAnalyticalSolution)
        {
          header << "#  (1)time        (2)dt             (3)#iterations   (4)[pot] L2 error (5)[pot] max error";
          for(int i=0; i<NUMBER_OF_SPECIES; ++i)
          {
            header << " (" << (count++) << ")[" << ION_NAMES[i] << "] L2 error";
            header << " (" << (count++) << ")[" << ION_NAMES[i] << "] max error";
          }
        } else {
          header << "#  (1)time        (2)dt             (3)#iterations";
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
              it != debugData.end(); ++it)
        {
          header << "   (" << (count++) << ")" << it->first;
        }
        return header.str();
      }

      void registerDebugData(std::string key, Real value = std::numeric_limits<Real>::lowest(),
          const int reductionFunction = ReductionFunction::Function_Identity)
      {
        debugData[key] = value;
        debugDataReductionFunctions[key] = reductionFunction;
      }

      void setDebugData(std::string key, Real value)
      {
        if(debugData.count(key))

          debugData[key] = value;
        else
          DUNE_THROW(Dune::Exception, "Key '"  << key << "' is not present in DiagnosticInfo::debugData map!");
      }

      std::vector<Real> getVectorRepresentation() const
      {
        int count = 0;

        std::vector<Real> vec;
        //vec[0] = time;
        vec.push_back(dt);
        vec.push_back((Real) iterations);
        if(hasAnalyticalSolution)
        {
          vec.push_back(l2ErrorPot);
          vec.push_back(maxErrorPot);
          for(int i=0; i<l2ErrorCon.size(); i++)
          {
            vec.push_back(l2ErrorCon[i]);
            vec.push_back(maxErrorCon[i]);
          }
        }
        for(typename std::map<std::string,Real>::const_iterator it = debugData.begin();
            it != debugData.end(); ++it)
        {
          vec.push_back(it->second);
        }
        return vec;
      }

      //! This returns a vector of the same size as getVectorRepresentation containing a flag
      //! for each diagnostic data entry which function to be used when reducing the values
      //! on all processors to the root node.
      std::vector<int> getReductionFunctions() const
      {
        int count = 0;

        std::vector<int> vec;
        //vec[0] = time;
        vec.push_back(ReductionFunction::Function_Identity); // dt
        vec.push_back(ReductionFunction::Function_Identity); // iterations
        if(hasAnalyticalSolution)
        {
          vec.push_back(ReductionFunction::Function_Identity); // l2ErrorPot
          vec.push_back(ReductionFunction::Function_Max); // maxErrorPot
          for(int i=0; i<l2ErrorCon.size(); i++)
          {
            vec.push_back(ReductionFunction::Function_Identity); // l2ErrorCon[i]
            vec.push_back(ReductionFunction::Function_Max); // maxErrorCon[i]
          }
        }
        for(typename std::map<std::string,int>::const_iterator it = debugDataReductionFunctions.begin();
            it != debugDataReductionFunctions.end(); ++it)
        {
          vec.push_back(it->second);
        }
        return vec;

      }

      void clear()
      {
        iterations = 0;
        dt = 0.0;
        l2ErrorPot = 0.0;
        maxErrorPot = 0.0;
        l2ErrorCon.clear(); l2ErrorCon.resize(NUMBER_OF_SPECIES);
        maxErrorCon.clear(); maxErrorCon.resize(NUMBER_OF_SPECIES);

        maxDiffCon.clear(); maxDiffCon.resize(0);
        maxDiffPot.clear(); maxDiffPot.resize(0);

        for(typename std::map<std::string, Real>::iterator it = debugData.begin();
            it != debugData.end(); ++it)
        {
          it->second = std::numeric_limits<Real>::lowest();
        }
      }
    };

    Acme2CylOutput(const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_,
        MultiGFS& multigfs_, GFS_CON& gfsCon_, GFS_POT& gfsPot_,
        Acme2CylSolutionVectors<Acme2CylTraits>& solutionVectors_,
        //U_CON& uCon_, U_POT& uPot_,
        typename Traits::BGF_MEMB_FLUX& bgfMembraneFlux_,
        typename Traits::BGF_MEMB_FLUX_MORI& bgfMoriFlux_,
        PHYSICS& physics_,
        typename Acme2CylTraits::PARAMETERS_POT& paramPot_,
        typename Acme2CylTraits::PARAMETERS_CON& paramCon_,
        int intorderPot_=2, int intorderCon_=2) :
      SUBSAMPLING_POINTS(0),
      vtkOutputBase_Domain("acme2_cyl"),
      vtkOutputBase_Membrane("acme2_cyl_membrane"),
      gv(gv_),
      elecGV(elecGV_),
      membGV(membGV_),
      multigfs(multigfs_),
      gfsCon(gfsCon_),
      gfsPot(gfsPot_),
      solutionVectors(solutionVectors_),
      u(solutionVectors_.unew),
      //uCon(solutionVectors_.unewCon),
      //uPot(solutionVectors_.unewPot),
      bgfMembraneFlux(bgfMembraneFlux_), // Reference (stateful object)!!
      bgfMoriFlux(bgfMoriFlux_), // Reference (stateful object)!!
      physics(physics_),
      intorderPot(intorderPot_),
      intorderCon(intorderCon_),
      totalInitParticles(NUMBER_OF_SPECIES),
      totalParticles(NUMBER_OF_SPECIES),
      sumIonFluxIntegrals(0.0),
      vtkwriterGV(gv,physics_.getParams().vtkSubSamplingLevel()),
      vtkwriterSubGV(membGV,0), // Don't do subsampling on the membrane!
      vtkBasenameDomain(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Domain),
      vtkBasenameMembrane(physics.getParams().getOutputPrefix() + "paraview/" + vtkOutputBase_Membrane),
      checkpointBasename(physics.getParams().getOutputPrefix() + "checkpoints/checkpoint"),
      hdf5Basename(physics.getParams().getOutputPrefix() + "hdf5/output"),
      fn_vtk_GV(vtkBasenameDomain),
      fn_vtk_SubGV(vtkBasenameMembrane),
      fn_checkpoints(checkpointBasename,-1), // checkpoints counter is always 1 lower than the others
      fn_hdf5(hdf5Basename),
      filePrefix(physics.getParams().getOutputPrefix()),
      dgfConElec(gfsCon,u),
      dgfConGradElec(gfsCon,u),
      dgfCon(gv,dgfConElec,physics),
      dgfConGrad(gv,dgfConGradElec,physics),
      dgfPot(gfsPot,u),
      dgfPotGrad(gfsPot,u),
      dgfPotGradElec(elecGV,dgfPotGrad,physics),
      gfChargeDensity(dgfCon,dgfCon,physics),
      gfConductivity(paramPot_,dgfCon,physics),
      gfElecFlux(dgfConElec,dgfConGradElec,dgfPotGradElec,physics),
      gfMembranePotentialMD(dgfPot, physics),
      gfMembranePotential(membGV, gfMembranePotentialMD, physics),
      gfMembraneFlux(membGV, bgfMembraneFlux, physics),
      bgfMembraneFlux_DiffTerm(bgfMembraneFlux),
      gfMembraneFlux_DiffTerm(membGV, bgfMembraneFlux_DiffTerm, physics),
      bgfMembraneFlux_DriftTerm(bgfMembraneFlux),
      gfMembraneFlux_DriftTerm(membGV, bgfMembraneFlux_DriftTerm, physics),
      bgfMembraneFlux_Leak(bgfMembraneFlux),
      gfMembraneFlux_Leak(membGV, bgfMembraneFlux_Leak, physics),
      bgfMembraneFlux_VoltageGated(bgfMembraneFlux),
      gfMembraneFlux_VoltageGated(membGV, bgfMembraneFlux_VoltageGated, physics),
      gfMembraneFlux_Mori(membGV, bgfMoriFlux, physics),
      bgfMembraneCurrent(paramCon_),
      bgfParamTotalMembraneFlux(paramCon_),
      bgfParamMembraneFluxWithoutMori(paramCon_, true),
      bgfPotNeumann(paramPot_,bgfParamTotalMembraneFlux),
      gfMembraneCurrent(membGV, bgfMembraneCurrent, physics),
      gfIonFlux(gv,gfElecFlux,bgfMembraneFlux,physics),
      gfAnalyticalSolutionCon(gv,physics.getParams()),
      gfAnalyticalSolutionPot(gv,physics.getParams()),
      vtkOutput(physics.getParams().doVTKOutput()),
      doFullGnuplotOutput(physics.getParams().doFullGnuplotOutput()),
      doHDF5Output(physics.getParams().doHDF5Output()),
      doFullHDF5Output(physics.getParams().doFullHDF5Output()),
      hdf5FileInitialized(false),
      to_mV(physics),
      timer(false),
      diagInfo(physics.getParams().hasAnalyticalSolution()),
      nTimeSteps(-1),
      useMori(physics.getParams().useMori())
    {
      timer.start();

      // All standard file I/O is handled by rank 0 only!
      // Exceptions: HDF5 output (MPI-parallel) and boundary value output (also rank size-1 for right boundary)
      if(gv.comm().rank() == 0)
      {
        std::vector<std::string> outputDirs = {"paraview", "checkpoints", "config", "hdf5"};
        int existingVTK = 0;
        int existingCheckpoints = 0;
        int existingHDF5 = 0;

        // Set up list of gnuplot files
        gnuplotFiles.push_back("cd.dat");
        gnuplotFiles.push_back("conductivity.dat");
        gnuplotFiles.push_back("pot.dat");
        gnuplotFiles.push_back("pot_grad.dat");
        gnuplotFiles.push_back("memb_pot.dat");

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          gnuplotFiles.push_back(physics.getIonName(j)+".dat");
          gnuplotFiles.push_back(physics.getIonName(j)+"_grad.dat");
          gnuplotFiles.push_back("memb_flux_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_flux_diff_term_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_flux_drift_term_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_flux_leak_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_flux_voltage_gated_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_flux_mori_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("memb_current_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("flux_" + physics.getIonName(j)+".dat");
          gnuplotFiles.push_back("con_error_" + physics.getIonName(j) + ".dat");
          gnuplotFiles.push_back("j_" + physics.getIonName(j)+".dat");

          gnuplotFiles.push_back("boundary_" + physics.getIonName(j) + "_bottom.dat");
          gnuplotFiles.push_back("boundary_" + physics.getIonName(j) + "_top.dat");
          gnuplotFiles.push_back("boundary_" + physics.getIonName(j) + "_left.dat");
          gnuplotFiles.push_back("boundary_" + physics.getIonName(j) + "_right.dat");
          gnuplotFiles.push_back("boundary_j_" + physics.getIonName(j) + "_bottom.dat");
          gnuplotFiles.push_back("boundary_j_" + physics.getIonName(j) + "_top.dat");
          gnuplotFiles.push_back("boundary_j_" + physics.getIonName(j) + "_left.dat");
          gnuplotFiles.push_back("boundary_j_" + physics.getIonName(j) + "_right.dat");

        }

        gnuplotFiles.push_back("boundary_pot_bottom.dat");
        gnuplotFiles.push_back("boundary_pot_top.dat");
        gnuplotFiles.push_back("boundary_pot_left.dat");
        gnuplotFiles.push_back("boundary_pot_right.dat");
        gnuplotFiles.push_back("boundary_j_pot_bottom.dat");
        gnuplotFiles.push_back("boundary_j_pot_top.dat");
        gnuplotFiles.push_back("boundary_j_pot_left.dat");
        gnuplotFiles.push_back("boundary_j_pot_right.dat");

        for(int k=0; k<physics.getMembrane().getChannelSet().size(); k++)
        {
          std::stringstream chStr;
          chStr << std::setw(2) << std::setfill('0') << k;
          gnuplotFiles.push_back("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat");
        }

        // 'One time' output stuff (only done once in the beginning)
        oneTimeGnuplotFiles.push_back("area_2d.dat");
        oneTimeGnuplotFiles.push_back("permittivity.dat");
        oneTimeGnuplotFiles.push_back("partition.dat");
        oneTimeGnuplotFiles.push_back("memb_groups.dat");
        oneTimeGnuplotFiles.push_back("volume.dat");
        oneTimeGnuplotFiles.push_back("base_surface.dat");
        oneTimeGnuplotFiles.push_back("side_surface.dat");
        oneTimeGnuplotFiles.push_back("surface.dat");


        // Check if output directory exists and create subdirs
        if(access(physics.getParams().getOutputPrefix().c_str(),0) != 0)
        {
          DUNE_THROW(Dune::IOError, "The specified output directory '" +
              physics.getParams().getOutputPrefix() + "' does not exist. Please create it first!");
        } else {

          // Create subdirs
          for(int i=0; i<outputDirs.size(); i++)
          {
            // Create single subdir and check for success
            if(access((physics.getParams().getOutputPrefix() + outputDirs[i]).c_str(),0) != 0)
            {
              int status = mkdir((physics.getParams().getOutputPrefix() + outputDirs[i]).c_str(),
                S_IRWXU | S_IRWXG | S_IROTH); // r+w+x access for user and group, r for other

              if(status != 0)
                DUNE_THROW(Dune::IOError,
                    "Could not create directory '" + (physics.getParams().getOutputPrefix() + outputDirs[i])
                    + "'");
            }
          }
        }

        // Clean up output dir only if the data shall not be appended to existing data files
        if(not (physics.getParams().doLoadState() && physics.getParams().doContinueSimulation()))
        {
          // Remove old VTK files
          std::string command = "rm " + vtkBasenameDomain + "*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());

          // Remove old checkpoint files
          command = "rm " + checkpointBasename + "*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());

          // Remove old VTK files
          command = "rm " + hdf5Basename + "*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());

          // Remove old files from config directory
          command = "rm " + physics.getParams().getOutputPrefix() + "config/*";
          debug_jochen << "Removing old files with command '" << command << "'" << std::endl;
          unusedVariable = std::system(command.c_str());

          // Following: Initialize (i.e., create and write headers) all gnuplot output files
          for(std::vector<std::string>::const_iterator it = gnuplotFiles.begin(); it != gnuplotFiles.end(); ++it )
          {
            Output::gnuplotInitialize(*it, filePrefix );
          }

        } else {
          // Continue simulation => keep old output files and count

          DIR *dir;
          struct dirent *ent;

          // Determine the number of written output VTK files
          dir = opendir ((physics.getParams().getOutputPrefix() + "paraview").c_str());
          if (dir != NULL) {

            /* print all the files and directories within directory */
            while ((ent = readdir(dir)) != NULL)
            {
              std::string filename(ent->d_name);
              if(filename.find(vtkOutputBase_Membrane) != std::string::npos)
              {
                existingVTK++;
              }
            }
            closedir (dir);
          } else {
            DUNE_THROW(Dune::IOError,
              "Error reading directory '" + (physics.getParams().getOutputPrefix() + "paraview")
              + "'");
          }

          // Determine the number of written output checkpoint files
          dir = opendir ((physics.getParams().getOutputPrefix() + "checkpoints").c_str());
          if (dir != NULL)
          {

            /* print all the files and directories within directory */
            while ((ent = readdir(dir)) != NULL)
            {
              std::string filename(ent->d_name);
              // Count number of checkpoint files generated by rank 0 only, each processor creates
              // its own checkpoint file!
              if(filename.find("checkpoint") != std::string::npos && filename.find("_p0") != std::string::npos)
              {
                existingCheckpoints++;
              }
            }
            closedir (dir);
          } else {
            DUNE_THROW(Dune::IOError,
              "Error reading directory '" + (physics.getParams().getOutputPrefix() + "checkpoints")
              + "'");
          }

          // Determine the number of written output checkpoint files
          dir = opendir ((physics.getParams().getOutputPrefix() + "hdf5").c_str());
          if (dir != NULL) {

            /* print all the files and directories within directory */
            while ((ent = readdir(dir)) != NULL)
            {
              std::string filename(ent->d_name);
              if(filename.find("output") != std::string::npos)
              {
                existingHDF5++;
              }
            }
            closedir (dir);
          } else {
            DUNE_THROW(Dune::IOError,
              "Error reading directory '" + (physics.getParams().getOutputPrefix() + "hdf5")
              + "'");
          }

          std::string command;
          // Backup one-time gnuplot files, as these are not handled by SimulationLoader::truncateAndDeleteExistingFiles()
          for(std::vector<std::string>::const_iterator it = oneTimeGnuplotFiles.begin(); it != oneTimeGnuplotFiles.end(); ++it )
          {
            command += "mv " + filePrefix + *it + " " + filePrefix + *it + ".old; ";
          }
          debug_jochen << "Backing up one-time gnuplot files with command '" << command << "'" << std::endl;
          int unusedVariable = std::system(command.c_str());

          // Calculate number of written time steps as the maximum of all three variants
          nExistingOutputFiles = std::max(existingVTK-1, std::max(existingCheckpoints, existingHDF5-1));
        }

        // One-time gnuplot files are always initialized (i.e. truncated), no matter if simulation is continued or not
        for(std::vector<std::string>::const_iterator it = oneTimeGnuplotFiles.begin(); it != oneTimeGnuplotFiles.end(); ++it )
        {
          Output::gnuplotInitialize(*it, filePrefix );
        }

        // Copy config file to config dir
        std::string command = "cp --backup=t " + physics.getParams().getConfigFileName() + " "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Copying config file using command '" << command << "'" << std::endl;
        int unusedVariable = std::system(command.c_str());

        // Create file containing the command line parameters used to start the program
        command = "echo \"" + physics.getParams().getCommandLineParameters()
            + "\" > command_line.dat && date >> command_line.dat "
            + "&& mv --backup=t command_line.dat "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Saving command line parameters using command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Generate detailed git version control info for this module
        command = "bash git-info.bash && mv --backup=t git-info.dat "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Saving Git VC info using command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Generate version control info for all Dune modules
        command = "bash dune-info.bash && mv --backup=t dune-info.dat "
            + physics.getParams().getOutputPrefix() + "config";
        debug_jochen << "Saving Dune VC info using command '" << command << "'" << std::endl;
        unusedVariable = std::system(command.c_str());

        // Create a DGF file from the multidomain GV
        if(physics.getParams().createGridFile() || physics.getParams().doEquilibration())
        {
          debug_jochen << "Saving DGF file..." << std::endl;
          Dune::DGFWriter<GV> dgfWriter(gv);

          std::string filename = physics.getParams().getOutputPrefix()
              + "config/" + physics.getParams().getSaveGridFileName();
          dgfWriter.write(filename);

#if USE_GRID == 1
          debug_jochen << "Saving special YaspGrid DGF file..." << std::endl;
          // Write out 'INTERVAL'-style DGF file in case of Yasp base grid
          int index = filename.find(".dgf");
          assert(index != std::string::npos);

          std::string basegrid_filename(filename.substr(0,index));
          basegrid_filename += "_basegrid.dgf";

          std::ofstream out_basegrid(basegrid_filename);
          if(! out_basegrid.good())
            DUNE_THROW(Dune::Exception, "DGF file '" << basegrid_filename << "' could not be created!");

          out_basegrid << "DGF" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "INTERVAL" << std::endl;
          out_basegrid << "0 0" << std::endl;
          out_basegrid << physics.getParams().xMax() << " " << physics.getParams().yMax() << std::endl;
          out_basegrid << (physics.getParams().X().size() - 1) << " "
                       << (physics.getParams().Y().size() - 1) << std::endl;
          out_basegrid << "# INTERVAL" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "GRIDPARAMETER" << std::endl;
          out_basegrid << "name Yasp_BaseGrid" << std::endl;
          out_basegrid << "overlap " << gv.overlapSize(0) << std::endl;
          out_basegrid << "# GRIDPARAMETER" << std::endl;
          out_basegrid << std::endl;
          out_basegrid << "# DGF";
#endif
        }
      }

      if(physics.getParams().general.get("oneTimeGridfunctionOutput",true))
      {
        debug_jochen << "Writing one-time gridfunction output..." << std::endl;

        // Write out volume gf's (this is done only once in the beginning)
        typename Traits::GF_PLAIN2D_AREA area2DGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, area2DGf, SUBSAMPLING_POINTS,
            filePrefix + "area_2d.dat", time, "plain 2D element area", true);

        // Write permittivity (gnuplot debug output)
        typename Traits::GF_PERMITTIVITY permittivityGf(gv,physics);
        GnuplotTools2D<Traits>::writeGF(physics, permittivityGf, SUBSAMPLING_POINTS,
            filePrefix + "permittivity.dat", time, "permittivity", true);

        typename Traits::GF_PARTITION processorGf(gv);
        GnuplotTools2D<Traits>::writeGF(physics, processorGf, SUBSAMPLING_POINTS,
              filePrefix + "partition.dat", time, "partition", true);

        if(physics.getParams().useMembrane())
        {
          // RangeType is tuple!
          typename Traits::GF_MEMB_GROUPS membGroupGf(gv,physics);

          std::string membGroupFilename = physics.getParams().getOutputPrefix() + "config/membrane_groups.dat";
          std::string infoMembGroups = "Membrane elements (1,2) coordinates (3) membrane index (4) group index (5) group name "
               "(6) permittivity (7) processor";
          infoMembGroups += "\n# " + physics.getGroupInfo();
          GnuplotTools2D<Traits>::writeGF(physics, membGroupGf, SUBSAMPLING_POINTS,
                filePrefix + "memb_groups.dat", time, infoMembGroups, true);
        }

#if HAVE_HDF5
        if(doHDF5Output)
        {
          std::string filename(physics.getParams().getOutputPrefix() + "hdf5/geometry.h5");
          HDF5Tools<Traits>::initialize(filename, 0.0,physics);
          std::string datasetName("area_2d");
          HDF5Tools<Traits>::writeGF(physics,area2DGf,SUBSAMPLING_POINTS,
              filename,datasetName);
          datasetName = "permittivity";
          HDF5Tools<Traits>::writeGF(physics,permittivityGf,SUBSAMPLING_POINTS,
              filename,datasetName);
          // FIXME Deactivated partition output, what is going wrong here?
          //datasetName = "partition";
          //HDF5Tools<Traits>::writeGF(physics,processorGf,SUBSAMPLING_POINTS,
          //    filename,datasetName);
          if(physics.getParams().useMembrane())
          {
            // This will currently not work, as the gridfunction returns a tuple of different data types
  //          datasetName = "membrane_element_indices";
  //          HDF5Tools<Traits>::writeGF(physics,membElementIndexGf,SUBSAMPLING_POINTS,
  //              filename,datasetName);
          }
        }
#endif

        // Helper grid functions for cylinder volume/surface areas etc.
        if(USE_CYLINDER_COORDINATES)
        {
          typename Traits::GF_VOLUME volumeGf(gv);
          GnuplotTools2D<Traits>::writeGF(physics, volumeGf, SUBSAMPLING_POINTS,
              filePrefix + "volume.dat", time, "cylinder element volume", true);

          typename Traits::GF_BASE_SURFACE baseSurfaceGf(gv);
          GnuplotTools2D<Traits>::writeGF(physics, baseSurfaceGf, SUBSAMPLING_POINTS,
              filePrefix + "base_surface.dat", time, "cylinder element base (top+bottom) areas", true);

          typename Traits::GF_SIDE_SURFACE sideSurfaceGf(gv);
          GnuplotTools2D<Traits>::writeGF(physics, sideSurfaceGf, SUBSAMPLING_POINTS,
              filePrefix + "side_surface.dat", time, "cylinder element side surface area", true);

          typename Traits::GF_SURFACE surfaceGf(gv);
          GnuplotTools2D<Traits>::writeGF(physics, surfaceGf, SUBSAMPLING_POINTS,
              filePrefix + "surface.dat", time, "total cylinder element surface", true);

#if HAVE_HDF5
          if(doHDF5Output)
          {
            std::string filename(physics.getParams().getOutputPrefix() + "hdf5/geometry.h5");
            std::string datasetName = "volume";
            HDF5Tools<Traits>::writeGF(physics,volumeGf,SUBSAMPLING_POINTS,
                filename,datasetName);
            datasetName = "base_surface";
            HDF5Tools<Traits>::writeGF(physics,baseSurfaceGf,SUBSAMPLING_POINTS,
                filename,datasetName);
            datasetName = "side_surface";
            HDF5Tools<Traits>::writeGF(physics,sideSurfaceGf,SUBSAMPLING_POINTS,
                filename,datasetName);
            datasetName = "surface";
            HDF5Tools<Traits>::writeGF(physics,surfaceGf,SUBSAMPLING_POINTS,
                filename,datasetName);
          }
#endif
        }
        debug_jochen << "...done!" << std::endl;
      }
      timer.stop();
    }

    void writeStep(double time_)
    {
      double start = timer.elapsed();
      double previous = start;
      double elapsed = start;
      timer.start();

      nTimeSteps++;
      time = time_;
      infoStream.str("");
      infoStream << "time: " << std::setprecision(16) << time;
      diagInfo.time = time;

#if HAVE_HDF5
      if(doHDF5Output)
      {
        if(! hdf5FileInitialized)
        {
          debug_jochen << "Initializing HDF file '" << getHDF5FileName() << "'" << std::endl;
          HDF5Tools<Traits>::initialize(getHDF5FileName(),time,physics);
          hdf5FileInitialized = true;
        }
        debug_jochen << "Writing to HDF5 file '" << getHDF5FileName() << "'" << std::endl;
      }
#endif

      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Output preparation time: " << (elapsed-previous) << " s" << std::endl;

      writePotentialOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Potential output time: " << (elapsed-previous) << " s" << std::endl;

      // Only do membrane output if there is a membrane!
      if(physics.getParams().useMembrane())
      {
        writeMembraneOutput();
        previous = elapsed;
        elapsed = timer.elapsed();
        debug_info << "[Acme2CylOutput] Membrane output time: " << (elapsed-previous) << " s" << std::endl;
      }

      writeConcentrationOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Concentration output time: " << (elapsed-previous) << " s" << std::endl;

      writeBoundaryOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Boundary output time: " << (elapsed-previous) << " s" << std::endl;

      writeDiagnosticOutput();
      previous = elapsed;
      elapsed = timer.elapsed();
      debug_info << "[Acme2CylOutput] Diagnostic output time: " << (elapsed-previous) << " s" << std::endl;

      if(vtkOutput)
      {
        //std::remove(fn_vtk_GV.getName());
        vtkwriterGV.write(fn_vtk_GV.getName(),Dune::VTK::OutputType::ascii);
        vtkwriterGV.clear();
        vtkwriterSubGV.write(fn_vtk_SubGV.getName(),Dune::VTK::OutputType::ascii);
        vtkwriterSubGV.clear();
        fn_vtk_GV.increment();
        fn_vtk_SubGV.increment();
      }
      if(doHDF5Output)
      {
        fn_hdf5.increment();
        hdf5FileInitialized = false;
      }
      if(physics.getParams().doCheckpointing())
      {
        fn_checkpoints.increment();
      }

      elapsed = timer.stop() - start;
      debug_info << "[Acme2CylOutput] TOTAL output time: " << elapsed << " s" << std::endl;
    }

    DiagnosticInfo& getDiagInfo()
    {
      return diagInfo;
    }

    void initDiagInfoFile()
    {
      if(gv.comm().rank() == 0)
      {
        Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
          physics.getParams().getOutputPrefix(), diagInfo.getHeader());
      }
    }


    void printAllCoeffs()
    {
      Output::printSingleCoefficientVector(u, "u");
    }

    double getTotalOutputTime()
    {
      return timer.elapsed();
    }

    int getNTimeSteps() const
    {
      return nTimeSteps;
    }

    int getIntOrderPot() const
    {
      return intorderPot;
    }

    int getIntOrderCon() const
    {
      return intorderCon;
    }

    /**
     * Function to write out a complete set of solution vectors which can be used to continue a
     * simulation later or to set a defined initial state for a new simulation
     *
     * @param time
     * @param dt
     * @param filename
     */
    void saveState(double time, double dt, std::string& filename)
    {
      timer.start();

      // Empty filename is given => this is a checkpoint output!
      if(filename.size() == 0)
      {
        std::stringstream checkpoint_filename;
        checkpoint_filename << fn_checkpoints.getName()
          << "_p" << physics.gridView().comm().rank()
          << ".dat";
        filename = checkpoint_filename.str();
      } else {
        int index = filename.find(".dat");
        assert(index != std::string::npos);
        std::stringstream save_filename;
        save_filename << filename.substr(0,index)
            << "_p" << physics.gridView().comm().rank()
            << ".dat";
        filename = save_filename.str();
      }

      Ax1SimulationData<Real> simulationData(gv.size(0), time, dt, filename);
      simulationData.setTimeStep(nTimeSteps);
      simulationData.setMpiRank(gv.comm().rank());
      simulationData.setMpiSize(gv.comm().size());

      // NEW method: save separate grid vectors for x,y coordinates
      simulationData.addVector("x", physics.getParams().X());
      simulationData.addVector("y", physics.getParams().Y());

      std::vector<Real> stdVec =  physics.getMembrane().getChannelSet().serializeChannelStates();
      simulationData.addVector("channelStates", stdVec);

      // Fill simulation data with std::vectors containing data from solutionVectors
      solutionVectors.serialize(simulationData);

      Ax1SimulationState<Real> state;
      state.setData(simulationData);
      state.saveState();

      timer.stop();
    }

    std::map<std::string,std::string> loadState(double& time, double& dt)
    {
      timer.start();

      // This class contains several helper functions when loading previously saved simulation data
      SimulationLoader<Acme2CylTraits, Acme2CylOutputTraits> simLoader(physics, gv, elecGV, gfsPot, gfsCon, solutionVectors);

      bool extrapolateXData = physics.getParams().doExtrapolateXData();

      std::map<std::string,std::string> loadFileNames = simLoader.getLoadFileNames();

      if(loadFileNames.size() > 2)
        DUNE_THROW(Dune::Exception, "Loading from more than 2 simulation states not implemented yet!");
      if(loadFileNames.size() < 1)
        DUNE_THROW(Dune::Exception, "Did not find any 'loadFilename' in config file!");

      bool doAllCoordinatesMatch = true;

      // Store all states we possibly load from in a vector and perform checks an them
      std::vector<Ax1SimulationState<Real> > states(loadFileNames.size());
      int i=0;
      for(std::map<std::string,std::string>::iterator it = loadFileNames.begin(); it != loadFileNames.end(); ++it)
      {
        debug_jochen << "Processing loadFilename " << it->second << " with group name '" << it->first << "'" << std::endl;

        // Get consistent parallel filename
        simLoader.getParallelFilename(it->second);

        if(gv.comm().size() == 1 && physics.getParams().general.get("loadFromParallelRun",false))
        {
          int np = physics.getParams().general.get("parallelRunNumberProcessors",1);
          states[i].loadState(gv,it->second,np);

        } else {
          // Load the simulation state and store it in a vector
          states[i].loadState(it->second);
        }

        debug_jochen << "Ok, loaded state: nElems = " << states[i].getData().getNElements() << std::endl;

        states[i].getData().addMetaInfo("groupName", it->first);

        // Check coordinates of the underlying grid; is it compatible?
        bool doCoordinatesMatch = simLoader.checkCoordinates(states[i].getData());

        doAllCoordinatesMatch = doAllCoordinatesMatch && doCoordinatesMatch;

        assert(doAllCoordinatesMatch == doCoordinatesMatch);

        i++;
      }


      // Get simulation data from first entry in states vector
      Ax1SimulationData<Real>& simulationData = states[0].getData();

      if(extrapolateXData && !(doAllCoordinatesMatch) )
      {
        // Load from vector of simulation states!
        simLoader.extrapolateData(states);

      } else {
        // If data didn't need to be extrapolated, simply copy loaded solution vectors!
        solutionVectors.deserialize(simulationData);

        // Load channel states
        if(physics.getParams().doLoadChannelStates())
        {
          if(not simulationData.hasVector("channelStates"))
            DUNE_THROW(Dune::Exception, "No vector 'channelStates' found in saved file!");

          if(extrapolateXData && !(doAllCoordinatesMatch))
          {
            debug_warn << "Cannot load channel states when grids don't match. Trying to load now..." << std::endl;
          }

          // Now load that shit!
          std::vector<Real> stdVec = simulationData.getVector("channelStates").data;
          physics.getMembrane().getChannelSet().deserializeChannelStates(stdVec);

          bgfMembraneFlux.markChannelsInitialized(true);
        }
      }
      time = simulationData.getTime();
      dt = simulationData.getDt();

      // Special case 'continueSimulation': Keep all data generated until the timestep we are loading from!
      if(physics.getParams().doLoadState() && physics.getParams().doContinueSimulation())
      {
        int loadTimeStep = simulationData.getTimeStep();

        // Set the timestep counter to the loaded value
        nTimeSteps = loadTimeStep;

        debug_jochen << "Continue simulation from timestep #" << loadTimeStep << " (t=" << time
            << ")" << std::endl;

        // Calculate number of written time steps as the maximum of all three variants
        nExistingOutputFiles = std::max(loadTimeStep, nExistingOutputFiles);


        std::vector<std::string> gnuplotFilesToBeTruncated = gnuplotFiles;
        gnuplotFilesToBeTruncated.push_back(physics.getParams().getDiagnosticsFilename());
        // Truncate gnuplot files / delete spurious checkpoint and HDF5 files
        simLoader.truncateAndDeleteExistingFiles(loadTimeStep, time, nExistingOutputFiles,
            gnuplotFilesToBeTruncated, vtkBasenameDomain, vtkBasenameMembrane,
            checkpointBasename, hdf5Basename, filePrefix);

        // Set filename helper counters to the correct value; do this for all processors!
        for(int i=0; i<loadTimeStep+1; i++)
        {
          fn_vtk_GV.increment();
          fn_vtk_SubGV.increment();
          fn_checkpoints.increment();
          fn_hdf5.increment();
        }
      }

      // simulationData is not valid after going out of this scope anyway!
      simulationData.clear();

      timer.stop();

      return loadFileNames;
    }


    /**
     * Call this to initialize a new output file, which can be used
     * to store data written out by the writeGFToGnuplot() method
     *
     * @param name
     */
    void initGnuplotFile(const std::string& name, const std::string header = "")
    {
      std::stringstream filename;
      filename << name << ".dat";
      Output::gnuplotInitialize(filename.str(), filePrefix, header);
    }


    /**
     * This is a method which other classes may use to write out non-standard gridfunctions.
     * This is flexible, but maybe a little more expensive (must check if file already exists,
     * gridfunction and solution vectors are created each time, etc.)
     *
     * Make sure to call initFile() once at the beginning to initialize the output file!
     *
     * @param gf
     * @param name
     */
    template<typename GF>
    void writeGFToGnuplot(GF& gf, const std::string& name, const std::string headerSuffix = "",
        const bool forceWrite = false)
    {
      GnuplotTools2D<Traits>::writeGF(physics, gf, SUBSAMPLING_POINTS, filePrefix + name + ".dat",
          time, infoStream.str() + headerSuffix, forceWrite);

#if HAVE_HDF5
      if(doHDF5Output)
      {
        if(! hdf5FileInitialized)
        {
          // As this method is called before writeStep(), the HDF5 has to be created first
          HDF5Tools<Traits>::initialize(getHDF5FileName(), diagInfo.time, physics);
          hdf5FileInitialized = true;
        }

        HDF5Tools<Traits>::writeGF(physics,gf,SUBSAMPLING_POINTS,getHDF5FileName(),name);
      }
#endif
    }

  private:
    void writePotentialOutput()
    {
      if(vtkOutput)
      {
        vtkwriterGV.addCellData(
            new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_POT>(dgfPot,"pot"));
        vtkwriterGV.addCellData(
            new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_POT_GRAD>(dgfPotGrad,"pot_grad"));
      }
      if(doFullGnuplotOutput)
      {
        GnuplotTools2D<Traits>::writeGF(physics, dgfPot, SUBSAMPLING_POINTS, filePrefix + "pot.dat",
            time, infoStream.str(), false, to_mV);

        // Get analytical solution for potential
        if(physics.getParams().hasAnalyticalSolution())
        {
          gfAnalyticalSolutionPot.setTime(time);
          GnuplotTools2D<Traits>::writeGF(physics, gfAnalyticalSolutionPot, SUBSAMPLING_POINTS,
              filePrefix + "pot_analytical.dat", time, infoStream.str(), false, to_mV);

          diagInfo.l2ErrorPot = ErrorNorms::l2Norm(dgfPot, gfAnalyticalSolutionPot, intorderPot);
          diagInfo.maxErrorPot = ErrorNorms::maxNorm(dgfPot, gfAnalyticalSolutionPot, intorderPot);
          debug_info << "== L2 error for potential: " << diagInfo.l2ErrorPot << std::endl;
          debug_info << "== Max error for potential: " << diagInfo.maxErrorPot << std::endl;
        }
        // potential gradient
        GnuplotTools2D<Traits>::writeGF(physics, dgfPotGrad, SUBSAMPLING_POINTS, filePrefix + "pot_grad.dat",
            time, infoStream.str());

        if(useMori)
        {
          GnuplotTools2D<Traits>::writeGF(physics, gfConductivity, SUBSAMPLING_POINTS, filePrefix + "conductivity.dat",
            time, infoStream.str(), false);
        }
      }



#if HAVE_HDF5
      if(doHDF5Output)
      {
        std::string datasetName("pot");
        HDF5Tools<Traits>::writeGF(physics,dgfPot,SUBSAMPLING_POINTS,
            getHDF5FileName(),datasetName, to_mV);
        datasetName = "pot_grad";
        HDF5Tools<Traits>::writeGF(physics,dgfPotGrad,SUBSAMPLING_POINTS,
            getHDF5FileName(),datasetName); // Do not convert gradient to mV, as
                                            // Matlab evaluation scripts already do the conversion

        if(useMori)
        {
          datasetName = "conductivity";
          HDF5Tools<Traits>::writeGF(physics,gfConductivity,SUBSAMPLING_POINTS,
              getHDF5FileName(),datasetName);
        }
      }
#endif
    }

    void writeMembraneOutput()
    {
      if(physics.getParams().useMembrane())
      {
        if(vtkOutput)
        {
          DUNE_THROW(Dune::Exception, "VTK output for membrane gridfunctions not implemented!");
          //vtkwriterSubGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_MEMB_POT>
          //  (gfMembranePotential,"memb_pot"));
        }

        GnuplotTools2D<Traits>::writeGF(physics, gfMembranePotentialMD, SUBSAMPLING_POINTS,
            filePrefix + "memb_pot.dat", time, std::string(""), false, to_mV);

#if HAVE_HDF5
        if(doHDF5Output)
        {
          std::string datasetName("memb_pot");
          HDF5Tools<Traits>::writeGF(physics,gfMembranePotentialMD,SUBSAMPLING_POINTS,
              getHDF5FileName(),datasetName,to_mV);
        }
#endif

        // Membrane flux
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          //debug_jochen << ">>> Getting flux component " << (j*GV::dimension + i) << std::endl;
          //componentFlux.push_back(j*GV::dimension);     // x component
          //int componentFlux = (j*GV::dimension + 1); // y component

          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX> gfSingleMembraneFlux(
              new typename Traits::GF_SINGLE_MEMB_FLUX(bgfMembraneFlux, j));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_DIFF> gfSingleMembraneFlux_DiffTerm(
              new typename Traits::GF_SINGLE_MEMB_FLUX_DIFF(bgfMembraneFlux_DiffTerm, j));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT> gfSingleMembraneFlux_DriftTerm(
              new typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT(bgfMembraneFlux_DriftTerm, j));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_LEAK> gfSingleMembraneFlux_Leak(
              new typename Traits::GF_SINGLE_MEMB_FLUX_LEAK(bgfMembraneFlux_Leak, j));
          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED> gfSingleMembraneFlux_VoltageGated(
              new typename Traits::GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED(bgfMembraneFlux_VoltageGated, j));

          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_FLUX_MORI> gfSingleMembraneFlux_Mori(
              new typename Traits::GF_SINGLE_MEMB_FLUX_MORI(bgfMoriFlux, j));

          Dune::shared_ptr<typename Traits::GF_SINGLE_MEMB_CURRENT> gfSingleMembraneCurrent(
              new typename Traits::GF_SINGLE_MEMB_CURRENT(bgfMembraneCurrent, j));

          Dune::shared_ptr<typename Traits::GF_SINGLE_PARAM_MEMB_FLUX> gfSingleParamTotalMembraneFlux(
              new typename Traits::GF_SINGLE_PARAM_MEMB_FLUX(bgfParamTotalMembraneFlux, j));
          Dune::shared_ptr<typename Traits::GF_SINGLE_PARAM_MEMB_FLUX> gfSingleParamMembraneFluxWithoutMori(
              new typename Traits::GF_SINGLE_PARAM_MEMB_FLUX(bgfParamMembraneFluxWithoutMori, j));

          if(vtkOutput)
          {
            DUNE_THROW(Dune::Exception, "VTK output for membrane gridfunctions not implemented!");
//            vtkwriterSubGV.addCellData(
//                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX>
//                        (gfSingleMembraneFlux,"memb_flux_" + physics.getIonName(j)));
//            vtkwriterSubGV.addCellData(
//                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_DIFF>
//                        (gfSingleMembraneFlux_DiffTerm,"memb_flux_diff_term_" + physics.getIonName(j)));
//            vtkwriterSubGV.addCellData(
//                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_DRIFT>
//                        (gfSingleMembraneFlux_DriftTerm,"memb_flux_drift_term_" + physics.getIonName(j)));
//            vtkwriterSubGV.addCellData(
//                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_LEAK>
//                        (gfSingleMembraneFlux_Leak,"memb_flux_leak_" + physics.getIonName(j)));
//            vtkwriterSubGV.addCellData(
//                new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_MEMB_FLUX_VOLTAGE_GATED>
//                        (gfSingleMembraneFlux_VoltageGated,"memb_flux_voltage_gated_" + physics.getIonName(j)));
          }

          if(bgfMembraneFlux.active())
          {
            // Default case: memb flux GF calculates the flux, use this for output
            GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux, SUBSAMPLING_POINTS,
              filePrefix + "memb_flux_" + physics.getIonName(j) + ".dat", time);
          } else {
            // In case the membrane flux is loaded from an external file, we must use the parameter class' j()
            // function to evaluate the actual membrane flux via an adapter class. We need a different adapter
            // class for Mori, as we want to exclude the Mori flux contributions in that case.
            if(useMori)
            {
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleParamMembraneFluxWithoutMori, SUBSAMPLING_POINTS,
                filePrefix + "memb_flux_" + physics.getIonName(j) + ".dat", time);
            } else {
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleParamTotalMembraneFlux, SUBSAMPLING_POINTS,
                filePrefix + "memb_flux_" + physics.getIonName(j) + ".dat", time);
            }
          }

          // In any case, it is useful to write out the actual value used in the operator from the param class'
          // j() function, even if this means some overhead when not using Mori contributions.
          GnuplotTools2D<Traits>::writeGF(physics, *gfSingleParamTotalMembraneFlux, SUBSAMPLING_POINTS,
              filePrefix + "j_" + physics.getIonName(j) + ".dat", time);

          if(doFullGnuplotOutput)
          {
            if(bgfMembraneFlux.active())
            {
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux_DiffTerm, SUBSAMPLING_POINTS,
                  filePrefix + "memb_flux_diff_term_" + physics.getIonName(j) + ".dat", time);
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux_DriftTerm, SUBSAMPLING_POINTS,
                  filePrefix + "memb_flux_drift_term_" + physics.getIonName(j) + ".dat", time);
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux_Leak, SUBSAMPLING_POINTS,
                  filePrefix + "memb_flux_leak_" + physics.getIonName(j) + ".dat", time);
              GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux_VoltageGated, SUBSAMPLING_POINTS,
                  filePrefix + "memb_flux_voltage_gated_" + physics.getIonName(j) + ".dat", time);
            }
            GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneCurrent, SUBSAMPLING_POINTS,
                filePrefix + "memb_current_" + physics.getIonName(j) + ".dat", time);
          }
          if(useMori || physics.getParams().general.get("forceMoriFluxCalculation", false))
          {
            GnuplotTools2D<Traits>::writeGF(physics, *gfSingleMembraneFlux_Mori, SUBSAMPLING_POINTS,
                filePrefix + "memb_flux_mori_" + physics.getIonName(j) + ".dat", time);
          }

#if HAVE_HDF5
          if(doHDF5Output)
          {
            std::string datasetName = "memb_flux_" + physics.getIonName(j);
            if(bgfMembraneFlux.active())
            {
              HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux,SUBSAMPLING_POINTS,
                  getHDF5FileName(),datasetName);
            } else {
              if(useMori)
              {
                HDF5Tools<Traits>::writeGF(physics,*gfSingleParamMembraneFluxWithoutMori,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
              } else {
                HDF5Tools<Traits>::writeGF(physics,*gfSingleParamTotalMembraneFlux,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
              }
            }
            datasetName = "j_" + physics.getIonName(j);
            HDF5Tools<Traits>::writeGF(physics,*gfSingleParamTotalMembraneFlux,SUBSAMPLING_POINTS,
                getHDF5FileName(),datasetName);

            if(doFullHDF5Output)
            {
              if(bgfMembraneFlux.active())
              {
                datasetName = "memb_flux_diff_term_" + physics.getIonName(j);
                HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux_DiffTerm,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
                datasetName = "memb_flux_drift_term_" + physics.getIonName(j);
                HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux_DriftTerm,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
                datasetName = "memb_flux_leak_" + physics.getIonName(j);
                HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux_Leak,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
                datasetName = "memb_flux_voltage_gated_" + physics.getIonName(j);
                HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux_VoltageGated,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
              }
              if(physics.getParams().useMori() || physics.getParams().general.get("forceMoriFluxCalculation", false))
              {
                datasetName = "memb_flux_mori_" + physics.getIonName(j);
                HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneFlux_Mori,SUBSAMPLING_POINTS,
                    getHDF5FileName(),datasetName);
              }
              datasetName = "memb_current_" + physics.getIonName(j);
              HDF5Tools<Traits>::writeGF(physics,*gfSingleMembraneCurrent,SUBSAMPLING_POINTS,
                  getHDF5FileName(),datasetName);
            }
          }
#endif
        }

        // Channel stuff
        int nChannels = physics.getMembrane().getChannelSet().size();
        for(int k=0; k<nChannels; k++)
        {
          // Create grid function for the k-th channel
          Dune::shared_ptr<typename Traits::GF_CHANNEL> gfChannel(
              new typename Traits::GF_CHANNEL(gv,physics,k));
          if(vtkOutput)
          {
            DUNE_THROW(Dune::Exception, "VTK output for membrane gridfunctions not implemented!");
            //vtkwriterSubGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_CHANNEL>
            //          (gfChannel,"channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName()));
          }

          GnuplotTools2D<Traits>::writeGF(physics, *gfChannel, SUBSAMPLING_POINTS,
              filePrefix + "channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat", time);

#if HAVE_HDF5
          if(doHDF5Output && doFullHDF5Output)
          {
            std::string datasetName("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName());
            HDF5Tools<Traits>::writeGF(physics,*gfChannel,SUBSAMPLING_POINTS,
                getHDF5FileName(),datasetName);
          }
#endif
        }
      }
    }

    void writeConcentrationOutput()
    {
      if(doFullGnuplotOutput)
      {
        GnuplotTools2D<Traits>::writeGF(physics, gfChargeDensity, SUBSAMPLING_POINTS, filePrefix + "cd.dat",
            time, infoStream.str());
      }

      sumIonFluxIntegrals = 0.0;
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        // select 1 component from the concentrations vector DGF
        std::vector<size_t> component;
        component.push_back(j);

        std::vector<size_t> componentGradient;
        for(int i=0; i<GV::dimension; i++)
        {
          componentGradient.push_back(j*GV::dimension + i);
        }

        Dune::shared_ptr<typename Traits::DGF_SINGLE_CON_MD> dgfSingleCon(
            new typename Traits::DGF_SINGLE_CON_MD(dgfCon, component));
        Dune::shared_ptr<typename Traits::DGF_SINGLE_CON_GRAD_MD> dgfSingleConGrad(
            new typename Traits::DGF_SINGLE_CON_GRAD_MD(dgfConGrad, componentGradient));
        Dune::shared_ptr<typename Traits::GF_SINGLE_ION_FLUX> gfSingleIonFlux(
            new typename Traits::GF_SINGLE_ION_FLUX(gfIonFlux, componentGradient));

         // VTK output
        if(vtkOutput)
        {
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_SINGLE_CON_MD>
            (dgfSingleCon,physics.getIonName(j)));
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::DGF_SINGLE_CON_GRAD_MD>
            (dgfSingleConGrad,physics.getIonName(j) + "_grad"));
          vtkwriterGV.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<typename Traits::GF_SINGLE_ION_FLUX>
            (gfSingleIonFlux,"flux_" + physics.getIonName(j)));
        }

        if(doFullGnuplotOutput)
        {
          GnuplotTools2D<Traits>::writeGF(physics, *dgfSingleCon, SUBSAMPLING_POINTS,
              filePrefix + physics.getIonName(j) + ".dat", time, infoStream.str());

          // Get analytical solution for concentrations
          if(physics.getParams().hasAnalyticalSolution())
          {
            Dune::shared_ptr<typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON> gfSingleAnalyticalSolutionCon(
                new typename Traits::SINGLE_ANALYTICAL_SOLUTION_CON(gfAnalyticalSolutionCon, component));
            GnuplotTools2D<Traits>::writeGF(physics, gfAnalyticalSolutionCon, SUBSAMPLING_POINTS,
                filePrefix + physics.getIonName(j) + "_analytical.dat", time, infoStream.str());

            diagInfo.setl2ErrorCon(j, ErrorNorms::l2Norm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
                intorderCon));
            diagInfo.setMaxErrorCon(j, ErrorNorms::maxNorm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
                intorderCon));
            debug_info << "== L2 error for " << physics.getIonName(j) << " concentration: "
                << diagInfo.getl2ErrorCon(j) << std::endl;
            debug_info << "== Max error for " << physics.getIonName(j) << " concentration: "
                << diagInfo.getMaxErrorCon(j) << std::endl;
          }
          GnuplotTools2D<Traits>::writeGF(physics, *dgfSingleConGrad, SUBSAMPLING_POINTS,
              filePrefix + physics.getIonName(j)+"_grad.dat", time, infoStream.str());
          GnuplotTools2D<Traits>::writeGF(physics, *gfSingleIonFlux, SUBSAMPLING_POINTS,
              filePrefix + "flux_" + physics.getIonName(j)+".dat", time, infoStream.str());

          //totalParticles[j] = Tools::integral( x, con[j] );
          typename Traits::DGF_SINGLE_CON_MD::Traits::RangeType concIntegral;
          //const int intorderadd = (USE_CYLINDER_COORDINATES ? 2 : 0);
          const int intorderadd = 0;
          Acme2CylGeometryTools::integrateGridFunctionOverCylinder(*dgfSingleCon, concIntegral,
              intorderCon + intorderadd);
          if(time == 0)
          {
            totalInitParticles[j] = gv.comm().sum(concIntegral);
          } else {
            totalParticles[j] = gv.comm().sum(concIntegral);
          }

          if(gv.comm().rank() == 0 && physics.getParams().general.get("gnuplotOutput",true))
          {
            std::vector<Real> chargeChange = {totalParticles[j] / totalInitParticles[j] - 1.0};
            GnuplotTools2D<Traits>::gnuplotAddLine("con_error_" + physics.getIonName(j) + ".dat",
                                time, chargeChange, filePrefix);
          }
        }

        //typename GF_SINGLE_IONFLUX::Traits::RangeType ionFluxIntegral;
        //Acme2CylGeometryTools::integrateGridFunctionOverCylinder(*gfSingleIonFlux, ionFluxIntegral, intorderCon);
        //debug_verb << "Integral " << ION_NAMES[j] << " ion flux: " << ionFluxIntegral << std::endl;
        //sumIonFluxIntegrals += ionFluxIntegral;

#if HAVE_HDF5
        if(doHDF5Output && doFullHDF5Output)
        {
          std::string datasetName("flux_" + physics.getIonName(j));
          HDF5Tools<Traits>::writeGF(physics,*gfSingleIonFlux,SUBSAMPLING_POINTS,
              getHDF5FileName(),datasetName);
        }
#endif
      }

#if HAVE_HDF5
      if(doHDF5Output)
      {
        // new GFs using GridFunctionSubSpace instead of SelectComponentGridFunctionAdapter
        std::string filename = getHDF5FileName();
        // concentrations
        Dune::ForLoop<write_concentration,0,NUMBER_OF_SPECIES-1>::apply(gv,gfsCon,u,physics,
            SUBSAMPLING_POINTS,filename);
        // concentration gradients
        Dune::ForLoop<write_concentration_gradient,0,NUMBER_OF_SPECIES-1>::apply(gv,gfsCon,u,physics,
              SUBSAMPLING_POINTS,filename);

        // charge density
        std::string datasetName("cd");
        HDF5Tools<Traits>::writeGF(physics,gfChargeDensity,SUBSAMPLING_POINTS,
            getHDF5FileName(),datasetName);
      }
#endif

      // Legacy code: Print charge density integral (= total charge) for each subdomain separately
      if(doFullGnuplotOutput)
      {
        typename Traits::GF_CD::Traits::RangeType cdIntegral;
        const int intorderadd = (USE_CYLINDER_COORDINATES ? 2 : 0);
        Acme2CylGeometryTools::integrateGridFunctionOverCylinder(gfChargeDensity, cdIntegral,
            intorderCon + intorderadd);
        debug_verb << "^^^ TOTAL charge = " << cdIntegral << std::endl;

        Acme2CylGeometryTools::integrateGridFunctionOverCylinderSubdomain(gfChargeDensity, physics, CYTOSOL, cdIntegral,
            intorderCon + intorderadd);
        debug_info << "^^^^ [" << SUBDOMAIN_NAMES[CYTOSOL] << "] charge = " << cdIntegral << std::endl;
        Acme2CylGeometryTools::integrateGridFunctionOverCylinderSubdomain(gfChargeDensity, physics, ES, cdIntegral,
            intorderCon + intorderadd);
        debug_info << "^^^^ [" << SUBDOMAIN_NAMES[ES] << "] charge = " << cdIntegral << std::endl;
        Acme2CylGeometryTools::integrateGridFunctionOverCylinderSubdomain(gfChargeDensity, physics, MEMBRANE, cdIntegral,
            intorderCon + intorderadd);
        debug_info << "^^^^ [" << SUBDOMAIN_NAMES[MEMBRANE] << "] charge = " << cdIntegral << std::endl;
        if(std::abs(cdIntegral) > 0)
        {
          debug_warn << "WARNING: Charge density inside membrane is non-zero!" << std::endl;
        }
      }

      //debug_verb << "SUM integrals ion flux: " << sumIonFluxIntegrals << std::endl;
    }

    void writeBoundaryOutput()
    {
      bool useTimeDependentBoundaryValues = physics.getParams().boundary.get("useTimeDependentBoundaryValuesCon", false)
          || physics.getParams().boundary.get("useTimeDependentBoundaryValuesPot", false);

      writeBoundaryOutput(dgfPot, std::string("pot"));

      if(useTimeDependentBoundaryValues)
      {
        writeBoundaryOutput(bgfPotNeumann, std::string("j_pot"));
      }

      for(int j=0; j<NUMBER_OF_SPECIES; j++)
      {
        Dune::shared_ptr<typename Traits::DGF_SINGLE_CON_MD> dgfSingleCon(
          new typename Traits::DGF_SINGLE_CON_MD(dgfCon, j));
        std::string name = physics.getIonName(j);
        writeBoundaryOutput(*dgfSingleCon, name);

        if(useTimeDependentBoundaryValues)
        {
          Dune::shared_ptr<typename Traits::GF_SINGLE_PARAM_MEMB_FLUX> gfSingleParamTotalMembraneFlux(
            new typename Traits::GF_SINGLE_PARAM_MEMB_FLUX(bgfParamTotalMembraneFlux, j));
          name = "j_" + physics.getIonName(j);
          writeBoundaryOutput(*gfSingleParamTotalMembraneFlux, name);
        }
      }
    }

    template<typename GF>
    void writeBoundaryOutput(const GF& gf, std::string name)
    {
      if(doFullGnuplotOutput || physics.getParams().boundary.get("writeBoundaryOutput",false))
      {
        debug_info << "    [writeBoundaryOutput] Proccessing boundary GF '" << name  << "'..." << std::endl;

        std::vector<std::vector<typename PHYSICS::Element::Geometry::GlobalCoordinate> > pos_vec;
        std::vector<std::vector<typename GF::Traits::RangeType> > sol_vec;
        std::vector<typename Ax1Output2D<Traits>::SolutionVectorInfo> info;

        const int strategy= OutputStrategy<Traits,GF>::value;

        Ax1Output2D<Traits>::template evaluateGFOnBoundary<
          PHYSICS,GF,strategy>::
            getSolutionVector(physics,gf,pos_vec,sol_vec,info);
        //Ax1Output2D<Traits>::getBoundarySolutionVector(physics,gf,pos_vec,sol_vec,info);

        std::string prefix("");
        for(int boundary = 0; boundary < 4; boundary++)
        {
          std::string filename = filePrefix + "boundary_" + name + "_";
          switch(boundary)
          {
            case BoundaryPositions::BOUNDARY_BOTTOM :
            {
              filename += "bottom";
              break;
            }
            case BoundaryPositions::BOUNDARY_TOP :
            {
              filename += "top";
              break;
            }
            case BoundaryPositions::BOUNDARY_LEFT :
            {
              filename += "left";
              break;
            }
            case BoundaryPositions::BOUNDARY_RIGHT :
            {
              filename += "right";
              break;
            }
            default :
              DUNE_THROW(Dune::Exception, "Unexpected boundary position identifier!");
          }
          filename += ".dat";

          std::vector<typename Traits::DGF_POT::Traits::DomainType>& pos = pos_vec[boundary];
          std::vector<typename Traits::DGF_POT::Traits::RangeType>& sol = sol_vec[boundary];
          //debug_jochen << "My vector size: " << sol.size() << std::endl;

          const int rank = gv.comm().rank();
          const int size = gv.comm().size();

          //debug_jochen << "MPI rank: " << rank << std::endl;
          //debug_jochen << "MPI size: " << size << std::endl;

          if(boundary == BoundaryPositions::BOUNDARY_BOTTOM
              || boundary == BoundaryPositions::BOUNDARY_TOP)
          {
            std::vector<int> all_dimensions_x = physics.nElements_AllProcessors(0); // copy, may not be const in MPI_Gatherv
            std::vector<int> displacements = physics.nOffset_AllProcessors(0); // copy, may not be const in MPI_Gatherv
            int total_size_x = std::accumulate(all_dimensions_x.begin(),all_dimensions_x.end(),0);

            // Case node output: 1 additional element in vector
            if(strategy != OutputPoints::OUTPUT_MEMBRANE)
            {
              total_size_x += 1;
            }

            // Gather partial vectors on the root node
            std::vector<typename Traits::DGF_POT::Traits::DomainType> pos_all(rank == 0 ? total_size_x : 0);
            std::vector<typename Traits::DGF_POT::Traits::RangeType> sol_all(rank == 0 ? total_size_x : 0);

            // Parallel run: Communication is necessary!
            if(size > 1)
            {
              if(strategy == OutputPoints::OUTPUT_MEMBRANE)
              {
                // Intersection boundary output is done at element centers => partition is already non-overlapping.
                // New version which allows different numbers of elements in x-direction
                int status = MPI_Gatherv(&pos[0], info[boundary].dimensions[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                    &pos_all[0], &all_dimensions_x[0], &displacements[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                    0, gv.comm());

                status = MPI_Gatherv(&sol[0], info[boundary].dimensions[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                    &sol_all[0], &all_dimensions_x[0], &displacements[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                    0, gv.comm());

              } else {
                // Let k_p be the number of elements on processor p
                // Take vertices (0,...,k_0) on root processor p=0
                // Take vertices (1,...,k_p) on each processor p >= 1
                // => get non-overlapping vertex partition for boundary values
                if(rank == 0)
                {
                  pos_all[0] = pos[0];
                  sol_all[0] = sol[0];
                }

                // Old version with same number of elements on all processors
  //              //debug_jochen << "Communicating vectors; my vector size: " << sol.size()
  //              //    << " global vector size: " << sol_all.size() << std::endl;
  //              gv.comm().gather(&pos[1],&pos_all[1],pos.size()-1,0);
  //              gv.comm().gather(&sol[1],&sol_all[1],pos.size()-1,0);

                // New version which allows different numbers of elements in x-direction
                int status = MPI_Gatherv(&pos[1], info[boundary].dimensions[0]-1,
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                    &pos_all[1], &all_dimensions_x[0], &displacements[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::DomainType>::getType(),
                    0, gv.comm());

                status = MPI_Gatherv(&sol[1], info[boundary].dimensions[0]-1,
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                    &sol_all[1], &all_dimensions_x[0], &displacements[0],
                    Dune::MPITraits<typename Traits::DGF_POT::Traits::RangeType>::getType(),
                    0, gv.comm());
              }

            } else {
              pos_all = pos;
              sol_all = sol;
            }

            // Do actual output
            if(gv.comm().rank() == 0)
            {
              GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos_all,sol_all,prefix,infoStream.str());
            }
          } else {
            // Complete boundary (left or right) is present on rank 0 and rank size-1.
            // => just output them directly without communication
            if(boundary == BoundaryPositions::BOUNDARY_LEFT)
            {
              if(rank == 0)
              {
                assert(sol.size() > 0);
                GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos,sol,prefix,infoStream.str(),false);
              } else {
                assert(sol.size() == 0);
              }
            }
            if(boundary == BoundaryPositions::BOUNDARY_RIGHT)
            {
              if(rank == size-1)
              {
                assert(sol.size() > 0);
                GnuplotTools2D<Traits>::gnuplotAddBlock(filename,pos,sol,prefix,infoStream.str(),false);
              } else {
                assert(sol.size() == 0);
              }
            }
          }
        }

      }
    }

    void writeDiagnosticOutput()
    {
      // Diagnostic output
      std::vector<Real> vec = diagInfo.getVectorRepresentation();
      std::vector<int> reductionFuncs = diagInfo.getReductionFunctions();

      // Reduce diagnostic data from all processors on root node
      for(int i=0; i<vec.size(); i++)
      {
        switch(reductionFuncs[i])
        {
          case DiagnosticInfo::ReductionFunction::Function_Identity :
          {
            // Do nothing
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Min :
          {
            vec[i] = gv.comm().min(vec[i]);
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Max :
          {
            vec[i] = gv.comm().max(vec[i]);
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Mean :
          {
            vec[i] = gv.comm().sum(vec[i]);
            vec[i] /= gv.comm().size();
            break;
          }
          case DiagnosticInfo::ReductionFunction::Function_Sum :
          {
            vec[i] = gv.comm().sum(vec[i]);
            break;
          }
          default:
            DUNE_THROW(Dune::Exception, "Could not find a reduction function for diagnostic data"
                << " entry #" << i << "!");
        }
      }
      // Print reduced diagnostic data on root node
      if(gv.comm().rank() == 0)
      {
        //Output::printVector(vec);
        GnuplotTools2D<Traits>::gnuplotAddLine(physics.getParams().getDiagnosticsFilename(),
          diagInfo.time, vec, filePrefix);
      }
    }

#if HAVE_HDF5
    template <int i>
    struct write_concentration
    {
     template<typename GFS>
     static void apply(const GV& gv, const GFS &gfs, const U& u, const PHYSICS& physics,
         const int subsamplingPoints, const std::string filename)
     {
       typedef typename GFS_ION<i>::type SubGFS;
       SubGFS subGFS(gfs);
       typedef Dune::PDELab::DiscreteGridFunction<SubGFS,U> DGF_SUB_CON;
       DGF_SUB_CON dgfSubCon(subGFS,u);
       typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_SUB_CON, PHYSICS> DGF_SUB_CON_MD;
       DGF_SUB_CON_MD dgfSubConMD(gv,dgfSubCon,physics);

       std::string datasetName = physics.getIonName(i);
       HDF5Tools<Traits>::writeGF(physics,dgfSubConMD,subsamplingPoints,filename,datasetName);

       //debug_jochen << "Wrote GF #" << i << std::endl;
     }
    };

    template <int i>
    struct write_concentration_gradient
    {
     template<typename GFS>
     static void apply(const GV& gv, const GFS &gfs, const U& u, const PHYSICS& physics,
         const int subsamplingPoints, const std::string filename)
     {
       typedef typename GFS_ION<i>::type SubGFS;
       SubGFS subGFS(gfs);
       typedef Dune::PDELab::DiscreteGridFunctionGradient<SubGFS,U> DGF_SUB_CON_GRAD;
       DGF_SUB_CON_GRAD dgfSubConGrad(subGFS,u);
       typedef Ax1MultiDomainGridFunctionExtension<GV,DGF_SUB_CON_GRAD,PHYSICS> DGF_SUB_CON_GRAD_MD;
       DGF_SUB_CON_GRAD_MD dgfSubConGradMD(gv,dgfSubConGrad,physics);

       std::string datasetName(physics.getIonName(i) + "_grad");
       HDF5Tools<Traits>::writeGF(physics,dgfSubConGradMD,subsamplingPoints,filename,datasetName);
     }
    };

 #endif

    std::string getHDF5FileName()
    {
      std::stringstream filename;
      filename << fn_hdf5.getName() << ".h5";
      return filename.str();
    }

    // Constructor parameters
    const GV& gv;
    const SubGV& elecGV;
    const SubGV& membGV;

    MultiGFS& multigfs;
    GFS_CON& gfsCon;
    GFS_POT& gfsPot;

    Acme2CylSolutionVectors<Acme2CylTraits>& solutionVectors;
    U&       u;
    //U_CON&   uCon;
    //U_POT&   uPot;
    typename Traits::BGF_MEMB_FLUX&   bgfMembraneFlux;
    typename Traits::BGF_MEMB_FLUX_MORI&   bgfMoriFlux;
    PHYSICS& physics;
    const int intorderPot;
    const int intorderCon;

    std::vector<Real>                 totalInitParticles;
    std::vector<Real>                 totalParticles;
    Real                              sumIonFluxIntegrals;

    Dune::SubsamplingVTKWriter<GV>     vtkwriterGV;
    Dune::SubsamplingVTKWriter<SubGV>  vtkwriterSubGV;

    const std::string vtkBasenameDomain;
    const std::string vtkBasenameMembrane;
    const std::string checkpointBasename;
    const std::string hdf5Basename;
    Dune::PDELab::FilenameHelper      fn_vtk_GV;
    Dune::PDELab::FilenameHelper      fn_vtk_SubGV;
    Dune::PDELab::FilenameHelper      fn_checkpoints;
    Dune::PDELab::FilenameHelper      fn_hdf5;

    const std::string filePrefix;

    typename Traits::DGF_CON         dgfConElec;
    typename Traits::DGF_CON_GRAD    dgfConGradElec;
    typename Traits::DGF_CON_MD      dgfCon;
    typename Traits::DGF_CON_GRAD_MD dgfConGrad;

    typename Traits::DGF_POT           dgfPot;
    typename Traits::DGF_POT_GRAD      dgfPotGrad;
    typename Traits::DGF_POT_GRAD_ELEC dgfPotGradElec;

    typename Traits::GF_CD           gfChargeDensity;
    typename Traits::GF_CONDUCTIVITY gfConductivity;
    typename Traits::GF_ELEC_FLUX    gfElecFlux;

    typename Traits::GF_MEMB_POT_MD  gfMembranePotentialMD;
    typename Traits::GF_MEMB_POT     gfMembranePotential;

    typename Traits::GF_MEMB_FLUX       gfMembraneFlux;
    typename Traits::BGF_MEMB_FLUX_DIFF  bgfMembraneFlux_DiffTerm;
    typename Traits::GF_MEMB_FLUX_DIFF  gfMembraneFlux_DiffTerm;
    typename Traits::BGF_MEMB_FLUX_DRIFT bgfMembraneFlux_DriftTerm;
    typename Traits::GF_MEMB_FLUX_DRIFT gfMembraneFlux_DriftTerm;
    typename Traits::BGF_MEMB_FLUX_LEAK bgfMembraneFlux_Leak;
    typename Traits::GF_MEMB_FLUX_LEAK gfMembraneFlux_Leak;
    typename Traits::BGF_MEMB_FLUX_VOLTAGE_GATED bgfMembraneFlux_VoltageGated;
    typename Traits::GF_MEMB_FLUX_VOLTAGE_GATED gfMembraneFlux_VoltageGated;

    typename Traits::GF_MEMB_FLUX_MORI gfMembraneFlux_Mori;

    typename Traits::BGF_MEMB_CURRENT bgfMembraneCurrent;
    typename Traits::BGF_PARAM_MEMB_FLUX bgfParamTotalMembraneFlux;
    typename Traits::BGF_PARAM_MEMB_FLUX bgfParamMembraneFluxWithoutMori;

    typename Traits::BGF_PARAM_POT_NEUMANN bgfPotNeumann;

    typename Traits::GF_MEMB_CURRENT gfMembraneCurrent;

    typename Traits::GF_ION_FLUX     gfIonFlux;

    typename Traits::ANALYTICAL_SOLUTION_CON    gfAnalyticalSolutionCon;
    typename Traits::ANALYTICAL_SOLUTION_POT    gfAnalyticalSolutionPot;

    const bool vtkOutput;
    const bool doFullGnuplotOutput;
    const bool doHDF5Output;
    const bool doFullHDF5Output;

    bool hdf5FileInitialized;

    // Functor for converting pot to mV
    Ax1OutputFunctors::to_millivolt<PHYSICS> to_mV;

    Dune::Timer timer;
    DiagnosticInfo diagInfo;

    int nTimeSteps;
    double time;
    std::stringstream infoStream;

    int nExistingOutputFiles;
    std::vector<std::string> gnuplotFiles;
    std::vector<std::string> oneTimeGnuplotFiles;

    const bool useMori;

};



#endif /* DUNE_AX1_ACME2CYL_OUTPUT_HH */

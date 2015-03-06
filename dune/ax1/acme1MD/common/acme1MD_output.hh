/*
 * acme1MD_output.hh
 *
 *  Created on: Aug 11, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1MD_OUTPUT_HH
#define DUNE_AX1_ACME1MD_OUTPUT_HH

#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/function/selectcomponent.hh>

#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/error_norms.hh>
#include <dune/ax1/common/ionfluxgridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_boundaryfunction_membranefunction_adapter.hh>
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/ax1_multidomaingridfunction.hh>
#include <dune/ax1/common/ax1_multidomaingridfunctionextension.hh>
#include <dune/ax1/common/ax1_multidomaingridfunctionrestriction.hh>
#include <dune/ax1/common/channelgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunction.hh>
#include <dune/ax1/common/membranefluxgridfunctionhelper.hh>
#include <dune/ax1/common/membranepotentialgridfunction.hh>



template<class GV, class SubGV, class GFS_CON, class GFS_POT, class U_CON, class U_POT, class PHYSICS>
class Acme1MDOutput
{
  public:

    typedef typename GV::Grid::ctype Coord;
    typedef double Real;

    static const int SUBSAMPLING_POINTS = 2;
    static const bool vtkOutput = false;

    // DGFs
    typedef Dune::PDELab::DiscreteGridFunction        <GFS_POT,U_POT> DGF_POT;
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;
    typedef Ax1MultiDomainGridFunctionRestriction<SubGV,DGF_POT_GRAD,PHYSICS> DGF_POT_GRAD_ELEC;

    typedef Dune::PDELab::VectorDiscreteGridFunction  <GFS_CON,U_CON> DGF_CON;
    typedef Dune::PDELab::VectorDiscreteGridFunctionGradient<GFS_CON,U_CON> DGF_CON_GRAD;

    typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON, PHYSICS> DGF_CON_MD;
    typedef Ax1MultiDomainGridFunctionExtension<GV, DGF_CON_GRAD, PHYSICS> DGF_CON_GRAD_MD;

    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_MD,1> DGF_SINGLE_CON_MD;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_GRAD_MD,1> DGF_SINGLE_CON_GRAD_MD;

    typedef ChargeDensityGridFunction<DGF_CON_MD,PHYSICS> GF_CD;

    // Meta grud function for the electrolyte ion flux taking only gridfunctions defined on the elec subdomain
    typedef IonFluxGridFunction<DGF_CON, DGF_CON_GRAD, DGF_POT_GRAD_ELEC, PHYSICS> GF_ELEC_FLUX;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ELEC_FLUX,1> GF_SINGLE_ELEC_FLUX;

    typedef MembraneFluxGridFunction<DGF_CON_MD,DGF_POT,PHYSICS> BGF_MEMB_FLUX;
    typedef MembraneFluxGridFunction_DiffusionTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DIFF;
    typedef MembraneFluxGridFunction_DriftTerm<BGF_MEMB_FLUX> BGF_MEMB_FLUX_DRIFT;

    typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX,PHYSICS> GF_MEMB_FLUX;
    typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DIFF,PHYSICS> GF_MEMB_FLUX_DIFF;
    typedef Ax1BoundaryFunctionMembraneFunctionAdapter<SubGV,BGF_MEMB_FLUX_DRIFT,PHYSICS> GF_MEMB_FLUX_DRIFT;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX,1> GF_SINGLE_MEMB_FLUX;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX_DIFF,1> GF_SINGLE_MEMB_FLUX_DIFF;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_MEMB_FLUX_DRIFT,1> GF_SINGLE_MEMB_FLUX_DRIFT;

    typedef Ax1MultiDomainGridFunction<GV,PHYSICS,GF_ELEC_FLUX,GF_MEMB_FLUX> GF_ION_FLUX;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_ION_FLUX,1> GF_SINGLE_ION_FLUX;

    typedef MembranePotentialGridFunction<DGF_POT, PHYSICS> GF_MEMB_POT_MD;
    typedef Ax1MultiDomainGridFunctionRestriction<SubGV,GF_MEMB_POT_MD,PHYSICS> GF_MEMB_POT;

    typedef ChannelGridFunction<SubGV,Real,PHYSICS> GF_CHANNEL;


    // Analytical solutions
    typedef typename PHYSICS::Traits::SOLUTION_CON ANALYTICAL_SOLUTION_CON;
    typedef typename PHYSICS::Traits::SOLUTION_POT ANALYTICAL_SOLUTION_POT;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<ANALYTICAL_SOLUTION_CON,1> SINGLE_ANALYTICAL_SOLUTION_CON;


    struct DiagnosticInfo
    {
      DiagnosticInfo()
      : time(0.0),
        iterations(0), dt(0.0),
        l2ErrorPot(0.0), maxErrorPot(0.0),
        l2ErrorCon(NUMBER_OF_SPECIES), maxErrorCon(NUMBER_OF_SPECIES),
        tEquilibrium(0.0),
        maxDiffCon(0),
        maxDiffPot(0)
      {}

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

      std::map<std::string, Real> debugData;


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

      std::string getHeader(PHYSICS& physics)
      {
        std::stringstream header;
        int count = 6;
        header << "#  (1)time        (2)dt             (3)#iterations   (4)[pot] L2 error (5)[pot] max error";
        for(int i=0; i<NUMBER_OF_SPECIES; ++i)
        {
          header << " (" << (count++) << ")[" << physics.getIonName(i) << "] L2 error";
          header << " (" << (count++) << ")[" << physics.getIonName(i) << "] max error";
        }
        for(std::map<std::string,Real>::const_iterator it = debugData.begin();
              it != debugData.end(); ++it)
        {
          header << "   (" << (count++) << ")" << it->first;
        }
        return header.str();
      }

      void addDebugData(std::string key, Real value)
      {
        debugData[key] = value;
      }

      std::vector<Real> getVectorRepresentation() const
      {
        int count = 0;

        std::vector<Real> vec;
        //vec[0] = time;
        vec.push_back(dt);
        vec.push_back((Real) iterations);
        vec.push_back(l2ErrorPot);
        vec.push_back(maxErrorPot);
        for(int i=0; i<l2ErrorCon.size()*2; i=i+2)
        {
          vec.push_back(l2ErrorCon[i]);
          vec.push_back(maxErrorCon[i]);
        }
        for(std::map<std::string,Real>::const_iterator it = debugData.begin();
            it != debugData.end(); ++it)
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
      }
    };

    Acme1MDOutput(const GV& gv_, const SubGV& elecGV_, const SubGV& membGV_,
        GFS_CON& gfsCon_, GFS_POT& gfsPot_,
        U_CON& uCon_, U_POT& uPot_,
        BGF_MEMB_FLUX& bgfMembraneFlux_,
        PHYSICS& physics_,
        int intorderPot_=2, int intorderCon_=2) :
      gv(gv_),
      elecGV(elecGV_),
      membGV(membGV_),
      gfsCon(gfsCon_),
      gfsPot(gfsPot_),
      uCon(uCon_),
      uPot(uPot_),
      bgfMembraneFlux(bgfMembraneFlux_), // Reference!
      physics(physics_),
      intorderPot(intorderPot_),
      intorderCon(intorderCon_),
      con(NUMBER_OF_SPECIES),
      conGrad(NUMBER_OF_SPECIES),
      membraneFlux(NUMBER_OF_SPECIES),
      membraneFlux_DiffTerm(NUMBER_OF_SPECIES),
      membraneFlux_DriftTerm(NUMBER_OF_SPECIES),
      ionFlux(NUMBER_OF_SPECIES),
      solutionCon(NUMBER_OF_SPECIES),
      totalInitParticles(NUMBER_OF_SPECIES),
      totalParticles(NUMBER_OF_SPECIES),
      sumIonFluxIntegrals(0.0),
      vtkwriter(gv,Dune::VTK::DataMode::conforming),
      fn_vtk("paraview/acme1MD_Pk"),
      filePrefix(physics.getParams().getOutputPrefix()),
      dgfConElec(gfsCon,uCon),
      dgfConGradElec(gfsCon,uCon),
      dgfCon(gv,dgfConElec,physics),
      dgfConGrad(gv,dgfConGradElec,physics),
      dgfPot(gfsPot,uPot),
      dgfPotGrad(gfsPot,uPot),
      dgfPotGradElec(elecGV,dgfPotGrad,physics),
      gfChargeDensity(dgfCon,dgfCon,physics),
      gfElecFlux(dgfConElec,dgfConGradElec,dgfPotGradElec,physics),
      gfMembranePotentialMD(dgfPot, physics),
      gfMembranePotential(membGV, gfMembranePotentialMD, physics),
      gfMembraneFlux(membGV, bgfMembraneFlux, physics),
      bgfMembraneFlux_DiffTerm(bgfMembraneFlux),
      gfMembraneFlux_DiffTerm(membGV, bgfMembraneFlux_DiffTerm, physics),
      bgfMembraneFlux_DriftTerm(bgfMembraneFlux),
      gfMembraneFlux_DriftTerm(membGV, bgfMembraneFlux_DriftTerm, physics),
      gfIonFlux(gv,gfElecFlux,gfMembraneFlux,physics),
      gfAnalyticalSolutionCon(gv,physics.getParams()),
      gfAnalyticalSolutionPot(gv,physics.getParams())
    {
      /*
      typename DGF_POT::Traits::RangeType dgfPotRT;
      typename DGF_POT::Traits::RangeFieldType dgfPotRFT;
      typename DGF_POT_GRAD::Traits::RangeType dgfPotGradRT;
      typename DGF_POT_GRAD::Traits::RangeFieldType dgfPotGradRFT;
      //typename DGF_POT_GRAD::RF dgfPotGradRF;
      typename DGF_CON::Traits::RangeType dgfConRT;
      typename DGF_CON::Traits::RangeFieldType dgfConRFT;
      typename DGF_CON::RF dgfConRF;
      typename DGF_CON_GRAD::Traits::RangeType dgfConGradRT;
      typename DGF_CON_GRAD::Traits::RangeFieldType dgfConGradRFT;
      typename DGF_CON_GRAD::RF dgfConGradRF;

      debug_verb << "===== Acme1MDOutput TYPEINFO ====" << std::endl;
      debug_verb << "DGF_POT::Traits::RangeType = " << Tools::getTypeName(dgfPotRT) << std::endl;
      debug_verb << "DGF_POT::Traits::RangeFieldType = " << Tools::getTypeName(dgfPotRFT) << std::endl;
      debug_verb << "DGF_POT_GRAD::Traits::RangeType = " << Tools::getTypeName(dgfPotGradRT) << std::endl;
      debug_verb << "DGF_POT_GRAD::Traits::RangeFieldType = " << Tools::getTypeName(dgfPotGradRFT) << std::endl;
      //debug_verb << "DGF_POT_GRAD::RF = " << Tools::getTypeName(dgfPotGradRF) << std::endl<< std::endl;
      debug_verb << "DGF_CON::Traits::RangeType = " << Tools::getTypeName(dgfConRT) << std::endl;
      debug_verb << "DGF_CON::Traits::RangeFieldType = " << Tools::getTypeName(dgfConRFT) << std::endl;
      debug_verb << "DGF_CON::RF = " << Tools::getTypeName(dgfConRF) << std::endl<< std::endl;
      debug_verb << "DGF_CON_GRAD::Traits::RangeType = " << Tools::getTypeName(dgfConGradRT) << std::endl;
      debug_verb << "DGF_CON_GRAD::Traits::RangeFieldType = " << Tools::getTypeName(dgfConGradRFT) << std::endl;
      debug_verb << "DGF_CON_GRAD::RF = " << Tools::getTypeName(dgfConGradRF) << std::endl<< std::endl;


      typename GF_IONFLUX::Traits::RangeType v_kacknase;
      debug_jochen << "V_KACKNASE :" << Tools::getTypeName(v_kacknase) << std::endl;

      typename GF_SINGLE_IONFLUX::Traits::RangeType kacknase;
      debug_jochen << "KACKNASE :" << Tools::getTypeName(kacknase) << std::endl;
      */
    }

    void writeStep(double time_)
    {
      time = time_;
      infoStream.str("");
      infoStream << "time: " << time;
      diagInfo.time = time;

      if(time == 0)
      {
        Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
            physics.getParams().getOutputPrefix(), diagInfo.getHeader(physics));


        Output::gnuplotInitialize("debug_output.dat",
            physics.getParams().getOutputPrefix(), diagInfo.getHeader(physics));

        Output::gnuplotInitialize("operator_split_debug.dat",
            physics.getParams().getOutputPrefix());

        // gnuplot file initialization
        Output::gnuplotInitialize( "cd.dat", filePrefix );
        Output::gnuplotInitialize( "pot.dat", filePrefix);
        Output::gnuplotInitialize( "pot_grad.dat", filePrefix );
        Output::gnuplotInitialize( "memb_pot.dat", filePrefix );

        for(int j=0; j<NUMBER_OF_SPECIES; ++j)
        {
          // gnuplot
          Output::gnuplotInitialize( physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( physics.getIonName(j)+"_grad.dat", filePrefix);
          Output::gnuplotInitialize( "memb_flux_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "flux_" + physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( "con_error_" + physics.getIonName(j) + ".dat", filePrefix );
        }

        for(int k=0; k<physics.getMembrane().getChannelSet().size(); k++)
        {
          std::stringstream chStr;
          chStr << std::setw(2) << std::setfill('0') << k;
          Output::gnuplotInitialize("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat",
              filePrefix);
        }
      }

      typename GF_CD::Traits::RangeType cdIntegral;
      Dune::PDELab::integrateGridFunction(gfChargeDensity, cdIntegral, intorderCon);
      debug_verb << "^^^ TOTAL charge = " << cdIntegral << std::endl;

      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, CYTOSOL, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[CYTOSOL] << "] charge = " << cdIntegral << std::endl;
      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, ES, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[ES] << "] charge = " << cdIntegral << std::endl;
      Tools::integrateGridFunctionOverSubdomain(gfChargeDensity, physics, MEMBRANE, cdIntegral, intorderCon);
      debug_verb << "^^^^ [" << SUBDOMAIN_NAMES[MEMBRANE] << "] charge = " << cdIntegral << std::endl;
      if(std::abs(cdIntegral) > 0)
      {
        debug_warn << "WARNING: Charge density inside membrane is non-zero!" << std::endl;
      }

      writePotentialOutput();
      writeMembraneOutput();
      writeConcentrationOutput();


      // TODO **************************************************************************************
      std::vector<Real> vec = diagInfo.getVectorRepresentation();
      //debug_verb << "***" << vec[0] <<" - " << vec[1] << std::endl;
      Output::gnuplotMultiAppend(physics.getParams().getDiagnosticsFilename(), diagInfo.time, vec, filePrefix);

      // ============== WRITE EXACT POTENTIAL======================================================
      // nice hack!
      /*
      if (time == 9)
      {
        std::valarray<Real> source(position), potExact(position);
        source = chargeDensity*physics.getPoissonConstant();
        for (int i=0; i<position.size(); ++i)
        {
          potExact[i]=Tools::potentialExact(position, source, permittivity, position[i]);
        }
        //potExact = potExact - potExact.min();
        Output::gnuplotArray("potExact.dat", position, physics.convertTo_mV(potExact) );
      }
      */
      // =========================================================================================

      if(vtkOutput)
      {
        vtkwriter.write(fn_vtk.getName(),Dune::VTK::OutputType::ascii);
        vtkwriter.clear();
        fn_vtk.increment();
      }
    }

    DiagnosticInfo& getDiagInfo()
    {
      return diagInfo;
    }

    void printPotentialCoeffs()
    {
      Output::printSingleCoefficientVector(uPot, "uPot");
    }

    void printConcentrationCoeffs()
    {
      Output::printMultipleComponentCoefficientVector(uCon, NUMBER_OF_SPECIES);
    }

    void printAllCoeffs()
    {
      printConcentrationCoeffs();
      printPotentialCoeffs();
    }

    // This is of course not the right place to manipulate coefficient vectors;
    // it should be moved to Ax1Newton where it is actually needed
    template<typename GFS, typename U>
    void updateChildCoeffs(GFS& gfs, U& unew)
    {
      //TODO Update to new subgrid setup!

      //Tools::compositeToChildCoefficientVector(gfs, unew, uCon, 0);
      //Tools::compositeToChildCoefficientVector(gfs, unew, uPot, 1);
    }

  private:
    void writePotentialOutput()
    {
      if(vtkOutput)
      {
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF_POT>(dgfPot,"pot"));
      }

      // TODO getMultiGroupSolutionVector?
      Tools::getMultiGroupSolutionVector(physics, dgfPot, SUBSAMPLING_POINTS, x, pot);
      pot = physics.convertTo_mV(pot);

      // Get analytical solution for potential
      if(physics.getParams().hasAnalyticalSolution())
      {
        gfAnalyticalSolutionPot.setTime(time);
        gfAnalyticalSolutionCon.setTime(time);
        Tools::getSolutionVector(gfAnalyticalSolutionPot, SUBSAMPLING_POINTS, x, solutionPot);
        Output::gnuplotAppendDoubleArray("pot.dat", x, pot,
            physics.convertTo_mV(solutionPot), infoStream.str(), filePrefix);

        diagInfo.l2ErrorPot = ErrorNorms::l2Norm(dgfPot, gfAnalyticalSolutionPot, gv, intorderPot);
        diagInfo.maxErrorPot = ErrorNorms::maxNorm(dgfPot, gfAnalyticalSolutionPot, gv, intorderPot);

        debug_info << "== L2 error for potential: " << diagInfo.l2ErrorPot << std::endl;
        debug_info << "== Max error for potential: " << diagInfo.maxErrorPot << std::endl;


        /*
        debug_info << "== L2 error for concentrations vector: "
            << ErrorNorms::l2Norm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridview()) << std::endl;
        debug_info << "== Max error for concentrations vector: "
            << ErrorNorms::maxNorm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridview()) << std::endl;
        */
      } else {
        Output::gnuplotAppendArray("pot.dat", x, pot, infoStream.str(), filePrefix);
      }

      Tools::getMultiGroupSolutionVector(physics, dgfPotGrad, SUBSAMPLING_POINTS, x, potGrad);
      Output::gnuplotAppendArray("pot_grad.dat", x, potGrad, infoStream.str(), filePrefix);
    }

    void writeMembraneOutput()
    {
      if(physics.getParams().useMembrane())
      {
        // Membrane potential
        Tools::getMembraneSolutionVector(physics, gfMembranePotential, SUBSAMPLING_POINTS, x, membranePotential);
        assert(membranePotential.size() == 1); // In 1D, there is only one membrane element
        membranePotential = physics.convertTo_mV(membranePotential);

        for(int i=0; i<membranePotential.size(); i++)
        {
          Output::gnuplotAppend("memb_pot.dat", time, membranePotential[i], filePrefix);
        }

        // Membrane flux
        for(int j=0; j<NUMBER_OF_SPECIES; j++)
        {
          std::vector<size_t> component;
          component.push_back(j);

          debug_jochen << "FLUX:" << std::endl;
          Dune::shared_ptr<GF_SINGLE_MEMB_FLUX> gfSingleMembraneFlux(new GF_SINGLE_MEMB_FLUX(gfMembraneFlux, component));
          Tools::getMembraneSolutionVector(physics, *gfSingleMembraneFlux, SUBSAMPLING_POINTS, x, membraneFlux[j]);

          debug_jochen << "FLUX diff term:" << std::endl;
          Dune::shared_ptr<GF_SINGLE_MEMB_FLUX_DIFF> gfSingleMembraneFlux_DiffTerm(
              new GF_SINGLE_MEMB_FLUX_DIFF(gfMembraneFlux_DiffTerm, component));
          Tools::getMembraneSolutionVector(physics, *gfSingleMembraneFlux_DiffTerm, SUBSAMPLING_POINTS, x, membraneFlux_DiffTerm[j]);

          debug_jochen << "FLUX drift term:" << std::endl;
          Dune::shared_ptr<GF_SINGLE_MEMB_FLUX_DRIFT> gfSingleMembraneFlux_DriftTerm(
              new GF_SINGLE_MEMB_FLUX_DRIFT(gfMembraneFlux_DriftTerm, component));
          Tools::getMembraneSolutionVector(physics, *gfSingleMembraneFlux_DriftTerm, SUBSAMPLING_POINTS, x, membraneFlux_DriftTerm[j]);

          assert(membraneFlux[j].size() == 1); // In 1D, there is only one membrane element

          for(int i=0; i<membraneFlux[j].size(); i++)
          {
            std::vector<Real> fluxTerms = {membraneFlux[j][i], membraneFlux_DiffTerm[j][i], membraneFlux_DriftTerm[j][i]};
            // Print fluxes as well as diffusion/drift terms
            Output::gnuplotMultiAppend("memb_flux_" + physics.getIonName(j) + ".dat", time, fluxTerms, filePrefix);
            // Print flux only
            //Output::gnuplotAppend("memb_flux_" + physics.getIonName(j) + ".dat", time, membraneFlux[j][i], filePrefix);
          }
        }

        // Channel stuff
        int nChannels = physics.getMembrane().getChannelSet().size();
        for(int k=0; k<nChannels; k++)
        {
          // Create grid function for the k-th channel
          GF_CHANNEL gfChannel(membGV,physics,k);
          std::valarray<typename GF_CHANNEL::Traits::RangeType> channelInfo;
          Tools::getMembraneSolutionVector(physics, gfChannel, SUBSAMPLING_POINTS, x, channelInfo);

          std::vector<Real> flatChannelInfo;
          Tools::flattenVector(channelInfo, flatChannelInfo);
          //std::stringstream chStr;
          //chStr << std::setw(2) << std::setfill('0') << k;
          Output::gnuplotMultiAppend("channel_" + physics.getMembrane().getChannelSet().getChannel(k).getName() + ".dat",
              time, flatChannelInfo, filePrefix);
        }

      }
    }

    void writeConcentrationOutput()
    {
      Tools::getMultiGroupSolutionVector(physics, gfChargeDensity, SUBSAMPLING_POINTS, x, chargeDensity);//DG
      Output::gnuplotAppendArray("cd.dat", x, chargeDensity, infoStream.str(), filePrefix);//DG

      sumIonFluxIntegrals = 0.0;
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        // select 1 component from the concentrations vector DGF
        std::vector<size_t> component;
        component.push_back(j);

        Dune::shared_ptr<DGF_SINGLE_CON_MD>      dgfSingleCon    (new DGF_SINGLE_CON_MD(dgfCon, component));
        Dune::shared_ptr<DGF_SINGLE_CON_GRAD_MD> dgfSingleConGrad(new DGF_SINGLE_CON_GRAD_MD(dgfConGrad, component));
        Dune::shared_ptr<GF_SINGLE_ION_FLUX>     gfSingleIonFlux(new GF_SINGLE_ION_FLUX(gfIonFlux, component));

         // VTK output
        if(vtkOutput)
        {
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF_SINGLE_CON_MD>
            (dgfSingleCon,physics.getIonName(j)));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF_SINGLE_CON_GRAD_MD>
            (dgfSingleConGrad,"grad_" + physics.getIonName(j)));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<GF_SINGLE_ION_FLUX>
            (gfSingleIonFlux,"flux_" + physics.getIonName(j)));
        }

        Tools::getMultiGroupSolutionVector(physics, *dgfSingleCon, SUBSAMPLING_POINTS, x, con[j]);

        // Get analytical solution for concentrations
        if(physics.getParams().hasAnalyticalSolution())
        {
          Dune::shared_ptr<SINGLE_ANALYTICAL_SOLUTION_CON> gfSingleAnalyticalSolutionCon
          (new SINGLE_ANALYTICAL_SOLUTION_CON(gfAnalyticalSolutionCon, component));
          Tools::getMultiGroupSolutionVector(physics, *gfSingleAnalyticalSolutionCon, SUBSAMPLING_POINTS, x, solutionCon[j]);

          diagInfo.setl2ErrorCon(j, ErrorNorms::l2Norm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
              gv, intorderCon));
          diagInfo.setMaxErrorCon(j, ErrorNorms::maxNorm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
              gv, intorderCon));
          debug_info << "== L2 error for " << physics.getIonName(j) << " concentration: "
              << diagInfo.getl2ErrorCon(j) << std::endl;
          debug_info << "== Max error for " << physics.getIonName(j) << " concentration: "
              << diagInfo.getMaxErrorCon(j) << std::endl;
        }
        Tools::getMultiGroupSolutionVector(physics, *dgfSingleConGrad, SUBSAMPLING_POINTS, x, conGrad[j]);
        Tools::getMultiGroupSolutionVector(physics, *gfSingleIonFlux,  SUBSAMPLING_POINTS, x, ionFlux[j]);


        //if (con[j].min() < 0.0) DUNE_THROW( Dune::Exception, "Shit, " << ION_NAMES[j] << " concentration negative!" );

        // gnuplot
        if(physics.getParams().hasAnalyticalSolution())
        {
          Output::gnuplotAppendDoubleArray(physics.getIonName(j) + ".dat", x, con[j], solutionCon[j],
              infoStream.str(), filePrefix);
        } else {
          Output::gnuplotAppendArray(physics.getIonName(j)+".dat", x, con[j], infoStream.str(), filePrefix);
        }
        Output::gnuplotAppendArray(physics.getIonName(j)+"_grad.dat", x, conGrad[j], infoStream.str(), filePrefix);//DG
        Output::gnuplotAppendArray("flux_" + physics.getIonName(j)+".dat", x, ionFlux[j], infoStream.str(), filePrefix);//DG


        //totalParticles[j] = Tools::integral( x, con[j] );
        typename DGF_SINGLE_CON_MD::Traits::RangeType concIntegral;
        Dune::PDELab::integrateGridFunction(*dgfSingleCon, concIntegral, intorderCon);
        if(time == 0)
        {
          totalInitParticles[j] = concIntegral;
        } else {
          totalParticles[j] = concIntegral;
        }
        Output::gnuplotAppend("con_error_"+physics.getIonName(j)+".dat",
                              time, (totalParticles[j] / totalInitParticles[j] - 1.0), filePrefix);

        //typename GF_SINGLE_IONFLUX::Traits::RangeType ionFluxIntegral;
        //Dune::PDELab::integrateGridFunction(*gfSingleIonFlux, ionFluxIntegral, intorderCon);
        //debug_verb << "Integral " << ION_NAMES[j] << " ion flux: " << ionFluxIntegral << std::endl;
        //sumIonFluxIntegrals += ionFluxIntegral;
      }
      //debug_verb << "SUM integrals ion flux: " << sumIonFluxIntegrals << std::endl;
    }

    // Constructor parameters
    const GV& gv;
    const SubGV& elecGV;
    const SubGV& membGV;
    GFS_CON& gfsCon;
    GFS_POT& gfsPot;
    U_CON&   uCon;
    U_POT&   uPot;
    BGF_MEMB_FLUX&   bgfMembraneFlux;
    PHYSICS& physics;
    const int intorderPot;
    const int intorderCon;

    // (gnuplot) output arrays
    std::valarray<Real>               x;
    std::valarray<Real>               pot;
    std::valarray<Real>               potGrad;
    std::valarray<Real>               chargeDensity;
    std::valarray<Real>               membranePotential;
    std::vector<std::valarray<Real> > membraneFlux;
    std::vector<std::valarray<Real> > membraneFlux_DiffTerm;
    std::vector<std::valarray<Real> > membraneFlux_DriftTerm;
    std::vector<std::valarray<Real> > con;
    std::vector<std::valarray<Real> > conGrad;
    std::vector<std::valarray<Real> > ionFlux;

    std::vector<std::valarray<Real> > solutionCon;
    std::valarray<Real>               solutionPot;

    std::vector<Real>                 totalInitParticles;
    std::vector<Real>                 totalParticles;

    Real                              sumIonFluxIntegrals;

    Dune::VTKWriter<GV>               vtkwriter;
    Dune::PDELab::FilenameHelper      fn_vtk;

    const std::string filePrefix;

    DGF_CON         dgfConElec;
    DGF_CON_GRAD    dgfConGradElec;
    DGF_CON_MD      dgfCon;
    DGF_CON_GRAD_MD dgfConGrad;

    DGF_POT           dgfPot;
    DGF_POT_GRAD      dgfPotGrad;
    DGF_POT_GRAD_ELEC dgfPotGradElec;

    GF_CD           gfChargeDensity;
    GF_ELEC_FLUX    gfElecFlux;

    GF_MEMB_POT_MD  gfMembranePotentialMD;
    GF_MEMB_POT     gfMembranePotential;

    GF_MEMB_FLUX       gfMembraneFlux;
    BGF_MEMB_FLUX_DIFF  bgfMembraneFlux_DiffTerm;
    GF_MEMB_FLUX_DIFF  gfMembraneFlux_DiffTerm;
    BGF_MEMB_FLUX_DRIFT bgfMembraneFlux_DriftTerm;
    GF_MEMB_FLUX_DRIFT gfMembraneFlux_DriftTerm;

    GF_ION_FLUX     gfIonFlux;

    ANALYTICAL_SOLUTION_CON    gfAnalyticalSolutionCon;
    ANALYTICAL_SOLUTION_POT    gfAnalyticalSolutionPot;

    DiagnosticInfo diagInfo;

    double time;
    std::stringstream infoStream;
};

#endif /* DUNE_AX1_ACME1MD_OUTPUT_HH */

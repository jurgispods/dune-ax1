/*
 * acme1_output.hh
 *
 *  Created on: Aug 11, 2011
 *      Author: jpods
 */

#ifndef DUNE_AX1_ACME1_OUTPUT_HH
#define DUNE_AX1_ACME1_OUTPUT_HH

#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/function/selectcomponent.hh>

#include <dune/ax1/common/chargedensitygridfunction.hh>
#include <dune/ax1/common/constants.hh>
#include <dune/ax1/common/error_norms.hh>
#include <dune/ax1/common/ionfluxgridfunction.hh>
#include <dune/ax1/common/output.hh>
#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_discretegridfunction.hh>
#include <dune/ax1/common/membranepotentialgridfunction.hh>



template<class GFS_CON, class GFS_POT, class U_CON, class U_POT, class PHYSICS>
class Acme1Output
{
  public:

    typedef typename GFS_POT::Traits::GridViewType GV;
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;

    static const int SUBSAMPLING_POINTS = 0;
    static const bool vtkOutput = false;

    // DGFs
    typedef Dune::PDELab::DiscreteGridFunction        <GFS_POT,U_POT> DGF_POT;
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS_POT,U_POT> DGF_POT_GRAD;

    typedef Ax1VectorDiscreteGridFunction  <PHYSICS,GFS_CON,U_CON> DGF_CON;
    typedef Ax1VectorDiscreteGridFunctionGradient<PHYSICS,GFS_CON,U_CON> DGF_CON_GRAD;

    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON,1> DGF_SINGLE_CON;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<DGF_CON_GRAD,1> DGF_SINGLE_CON_GRAD;

    typedef ChargeDensityGridFunction<DGF_CON,PHYSICS> GF_CD;

    typedef IonFluxGridFunction<DGF_CON, DGF_CON_GRAD, DGF_POT_GRAD, PHYSICS> GF_IONFLUX;
    typedef Dune::PDELab::SelectComponentGridFunctionAdapter<GF_IONFLUX,1> GF_SINGLE_IONFLUX;

    typedef MembranePotentialGridFunction<DGF_POT, PHYSICS> GF_MEMB_POT;

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
        header << "#  time        dt             #iterations   [pot] L2 error [pot] max error";
        for(int i=0; i<NUMBER_OF_SPECIES; ++i)
        {
          header << " [" << physics.getIonName(i) << "] L2 error";
          header << " [" << physics.getIonName(i) << "] max error";
        }
        return header.str();
      }

      std::vector<Real> getVectorRepresentation() const
      {
        int size = 5;
        std::vector<Real> vec(size+l2ErrorCon.size()*2);
        vec[0] = time;
        vec[1] = dt;
        vec[2] = (Real) iterations;
        vec[3] = l2ErrorPot;
        vec[4] = maxErrorPot;
        for(int i=0; i<l2ErrorCon.size()*2; i=i+2)
        {
          vec[size+i] = l2ErrorCon[i];
          vec[size+i+1] = maxErrorCon[i];
        }
        return vec;
      }

      void registerDebugData(std::string key, Real value = -1e100)
      {
        DUNE_THROW(Dune::NotImplemented, "Acme0Output::DiagnosticInfo::registerDebugData");
      }

      void setDebugData(std::string key, Real value)
      {
        DUNE_THROW(Dune::NotImplemented, "Acme0Output::DiagnosticInfo::setDebugData");
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

    Acme1Output(GFS_CON& gfsCon_, GFS_POT& gfsPot_, U_CON& uCon_, U_POT& uPot_, PHYSICS& physics_,
        int intorderPot_=2, int intorderCon_=2) :
      gfsCon(gfsCon_),
      gfsPot(gfsPot_),
      physics(physics_),
      uCon(uCon_),
      uConPrevious(uCon_),
      uPot(uPot_),
      con(NUMBER_OF_SPECIES),
      conGrad(NUMBER_OF_SPECIES),
      ionFlux(NUMBER_OF_SPECIES),
      solutionCon(NUMBER_OF_SPECIES),
      totalInitParticles(NUMBER_OF_SPECIES),
      totalParticles(NUMBER_OF_SPECIES),
      sumIonFluxIntegrals(0.0),
      vtkwriter(gfsCon.gridView(),Dune::VTK::DataMode::conforming),
      fn_vtk("paraview/acme1_Pk"),
      filePrefix(physics.getParams().getOutputPrefix()),
      dgfCon(physics,gfsCon,uCon),
      dgfConPrevious(physics,gfsCon,uConPrevious),
      dgfConGrad(physics,gfsCon,uCon),
      dgfPot(gfsPot,uPot),
      dgfPotGrad(gfsPot,uPot),
      gfChargeDensity(dgfCon,dgfConPrevious,physics),
      gfIonFlux(dgfCon,dgfConGrad,dgfPotGrad,physics),
      gfMembranePotential(dgfPot, physics),
      gfAnalyticalSolutionCon(gfsCon.gridView(),physics.getParams()),
      gfAnalyticalSolutionPot(gfsPot.gridView(),physics.getParams()),
      intorderPot(intorderPot_),
      intorderCon(intorderCon_)
    {
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

      debug_verb << "===== Acme1Output TYPEINFO ====" << std::endl;
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

      typename GF_SINGLE_IONFLUX::Traits::RangeType kacknase;
      debug_jochen << "KACKNASE :" << Tools::getTypeName(kacknase) << std::endl;
    }

    void writeStep(double time)
    {
      bool useLogarithmicScaling = physics.getParams().useLogScaling();
      std::stringstream infoStream;
      infoStream << "time: " << time;
      diagInfo.time = time;

      if(time == 0)
      {
        Output::gnuplotInitialize(physics.getParams().getDiagnosticsFilename(),
            physics.getParams().getOutputPrefix(), diagInfo.getHeader(physics));

        Output::gnuplotInitialize("operator_split_debug.dat",
            physics.getParams().getOutputPrefix());

        // gnuplot file initialization
        Output::gnuplotInitialize( "cd.dat", filePrefix );
        Output::gnuplotInitialize( "pot.dat", filePrefix);
        Output::gnuplotInitialize( "pot_grad.dat", filePrefix );
        Output::gnuplotInitialize( "memb_pot.dat", filePrefix );
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

      // ============== POTENTIAL ===============================================================
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

        diagInfo.l2ErrorPot = ErrorNorms::l2Norm(dgfPot, gfAnalyticalSolutionPot, gfsPot.gridView(), intorderPot);
        diagInfo.maxErrorPot = ErrorNorms::maxNorm(dgfPot, gfAnalyticalSolutionPot, gfsPot.gridView(), intorderPot);

        debug_info << "== L2 error for potential: " << diagInfo.l2ErrorPot << std::endl;
        debug_info << "== Max error for potential: " << diagInfo.maxErrorPot << std::endl;


        /*
        debug_info << "== L2 error for concentrations vector: "
            << ErrorNorms::l2Norm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridView()) << std::endl;
        debug_info << "== Max error for concentrations vector: "
            << ErrorNorms::maxNorm(dgfCon, gfAnalyticalSolutionCon, gfsCon.gridView()) << std::endl;
        */
      } else {
        Output::gnuplotAppendArray("pot.dat", x, pot, infoStream.str(), filePrefix);
      }

#if USE_SUBGRID==1 || USE_SUBGRID==2
      // Get potential jump at membrane (1D only!)
      if(physics.getParams().useMembrane())
      {
        Tools::getMultiGroupSolutionVector(physics, gfMembranePotential, SUBSAMPLING_POINTS, x, membranePotential);
        membranePotential = physics.convertTo_mV(membranePotential);
        Output::gnuplotAppendArray("memb_pot.dat", x, membranePotential, infoStream.str(), filePrefix);
      }
#endif

      /*
      if(physics.getParams().useMembrane())
      {
        Real pot_l, pot_r;
        for(int i=0; i<x.size(); ++i)
        {
          if (std::abs(std::abs(x[i]) - 0.5*physics.getParams().dMemb()) < 1e-12)
          {
            if(x[i] < 0)
            {
              pot_l = pot[i];
              i++;
            }
            if(x[i] > 0)
            {
              pot_r = pot[i];
              break;
            }
          }
        }
        debug_info << " ++ Potential jump at membrane: " << (pot_r - pot_l) << std::endl;
      }
      */

      // TODO getMultiGroupSolutionVector?
      Tools::getMultiGroupSolutionVector(physics, dgfPotGrad, SUBSAMPLING_POINTS, x, potGrad);
      Output::gnuplotAppendArray("pot_grad.dat", x, potGrad, infoStream.str(), filePrefix);
      // =========================================================================================

      // TODO getMultiGroupSolutionVector?
      Tools::getMultiGroupSolutionVector(physics, gfChargeDensity, SUBSAMPLING_POINTS, x, chargeDensity);//DG
      Output::gnuplotAppendArray("cd.dat", x, chargeDensity, infoStream.str(), filePrefix);//DG

      sumIonFluxIntegrals = 0.0;
      // ============== CONCENTRATION ============================================================
      for(int j=0; j<NUMBER_OF_SPECIES; ++j)
      {
        if(time == 0)
        {
          // gnuplot
          Output::gnuplotInitialize( physics.getIonName(j)+".dat", filePrefix);
          Output::gnuplotInitialize( physics.getIonName(j)+"_grad.dat", filePrefix);
          Output::gnuplotInitialize( physics.getIonName(j)+"_flux.dat", filePrefix);

          Output::gnuplotInitialize( "con_error_" + physics.getIonName(j) + ".dat", filePrefix );
        }

        // select 1 component from the concentrations vector DGF
        std::vector<size_t> component;
        component.push_back(j);

        Dune::shared_ptr<DGF_SINGLE_CON>      dgfSingleCon    (new DGF_SINGLE_CON(dgfCon, component));
        Dune::shared_ptr<DGF_SINGLE_CON_GRAD> dgfSingleConGrad(new DGF_SINGLE_CON_GRAD(dgfConGrad, component));
        Dune::shared_ptr<GF_SINGLE_IONFLUX>   gfSingleIonFlux(new GF_SINGLE_IONFLUX(gfIonFlux, component));

         // VTK output
        if(vtkOutput)
        {
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF_SINGLE_CON>
            (dgfSingleCon,physics.getIonName(j)));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF_SINGLE_CON_GRAD>
            (dgfSingleConGrad,"grad_" + physics.getIonName(j)));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<GF_SINGLE_IONFLUX>
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
              gfsCon.gridView(), intorderCon));
          diagInfo.setMaxErrorCon(j, ErrorNorms::maxNorm(*dgfSingleCon, *gfSingleAnalyticalSolutionCon,
              gfsCon.gridView(), intorderCon));
          debug_info << "== L2 error for " << physics.getIonName(j) << " concentration: "
              << diagInfo.getl2ErrorCon(j) << std::endl;
          debug_info << "== Max error for " << physics.getIonName(j) << " concentration: "
              << diagInfo.getMaxErrorCon(j) << std::endl;
        }
        Tools::getMultiGroupSolutionVector(physics, *dgfSingleConGrad, SUBSAMPLING_POINTS, x, conGrad[j]);
        Tools::getMultiGroupSolutionVector(physics, *gfSingleIonFlux,  SUBSAMPLING_POINTS, x, ionFlux[j]);

        if(useLogarithmicScaling)
        {
          con[j] = std::exp(con[j]);
          conGrad[j] = std::exp(conGrad[j]);
          DUNE_THROW(Dune::NotImplemented, "Logarithmic scaling for ion flux not implemented!");
        }

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
        Output::gnuplotAppendArray(physics.getIonName(j)+"_flux.dat", x, ionFlux[j], infoStream.str(), filePrefix);//DG


        //totalParticles[j] = Tools::integral( x, con[j] );
        typename DGF_SINGLE_CON::Traits::RangeType concIntegral;
        Dune::PDELab::integrateGridFunction(*dgfSingleCon, concIntegral, intorderCon);
        if(time == 0)
        {
          totalInitParticles[j] = concIntegral;
        } else {
          totalParticles[j] = concIntegral;
        }
        Output::gnuplotAppend("con_error_"+physics.getIonName(j)+".dat",
                              time, (totalParticles[j] / totalInitParticles[j] - 1.0), filePrefix);

        typename GF_SINGLE_IONFLUX::Traits::RangeType ionFluxIntegral;
        Dune::PDELab::integrateGridFunction(*gfSingleIonFlux, ionFluxIntegral, intorderCon);
        //debug_verb << "Integral " << ION_NAMES[j] << " ion flux: " << ionFluxIntegral << std::endl;
        sumIonFluxIntegrals += ionFluxIntegral;
      }
      //debug_verb << "SUM integrals ion flux: " << sumIonFluxIntegrals << std::endl;
      // =========================================================================================

      // TODO **************************************************************************************
      std::vector<Real> vec = diagInfo.getVectorRepresentation();
      //debug_verb << "***" << vec[0] <<" - " << vec[1] << std::endl;
      Output::gnuplotMultiAppend(physics.getParams().getDiagnosticsFilename(), diagInfo.time,
        vec, "", filePrefix);

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
      uConPrevious = uCon;

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
      Output::printSingleCoefficientVectorDG(uPot, "pot");
    }

    void printConcentrationCoeffs()
    {
      Output::printMultipleComponentCoefficientVectorDG(uCon, NUMBER_OF_SPECIES);
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
    GFS_CON& gfsCon;
    GFS_POT& gfsPot;
    PHYSICS& physics;
    U_CON&   uCon;
    U_CON    uConPrevious;
    U_POT&   uPot;

    // (gnuplot) output arrays
    std::valarray<Real>               x;
    std::valarray<Real>               pot;
    std::valarray<Real>               potGrad;
    std::valarray<Real>               chargeDensity;
    std::valarray<Real>               membranePotential;
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

    DGF_CON         dgfCon;
    DGF_CON         dgfConPrevious;
    DGF_CON_GRAD    dgfConGrad;
    DGF_POT         dgfPot;
    DGF_POT_GRAD    dgfPotGrad;
    GF_CD           gfChargeDensity;
    GF_IONFLUX      gfIonFlux;
    GF_MEMB_POT     gfMembranePotential;

    ANALYTICAL_SOLUTION_CON    gfAnalyticalSolutionCon;
    ANALYTICAL_SOLUTION_POT    gfAnalyticalSolutionPot;

    DiagnosticInfo diagInfo;

    const int intorderPot;
    const int intorderCon;
};


#endif /* DUNE_AX1_ACME1_OUTPUT_HH */

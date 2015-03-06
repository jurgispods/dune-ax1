
#ifndef DUNE_AX1_CONSTANTS_HH
#define DUNE_AX1_CONSTANTS_HH

#include <string>
#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/debugstream.hh>

#include <dune/ax1/common/ax1_parallelhelper.hh>

#if USE_GRID==1
typedef Dune::YaspGrid<2,Dune::TensorProductCoordinates<double,2> > BaseGrid;
#elif USE_GRID==2
typedef Dune::UGGrid<2> BaseGrid;
#endif

#ifdef MULTIPLE_MEMBRANE_ELEMENTS
#define DO_MULTIDOMAIN_BLOCKING (0)
#else
#define DO_MULTIDOMAIN_BLOCKING (1)
#endif

//! \brief Number of ion species
static const unsigned int NUMBER_OF_SPECIES = 3;


//! \brief Blocksize in ISTL Backends
#if DO_MULTIDOMAIN_BLOCKING==1
// Use this for a fully blocked matrix
static const unsigned int AX1_BLOCKSIZE = NUMBER_OF_SPECIES+1;
#else
// Use this for a non-blocked matrix
static const unsigned int AX1_BLOCKSIZE = 1;
#endif

//! \brief Flag for using cylinder coordinates in acme2_cyl
#ifdef USE_CYLINDER_COORDS
static const bool USE_CYLINDER_COORDINATES = true;
#else
static const bool USE_CYLINDER_COORDINATES = false;
#endif

//! \brief enum for ion species
enum { Na = 0, K = 1, Cl = 2, Ca = 3 };
static const std::vector<std::string> ION_NAMES = { "Na", "K", "Cl"};

//! \brief enum for element subdomains (formerly known as 'element groups')
enum { CYTOSOL = 0, ES = 1, MEMBRANE = 2};
static const std::vector<std::string> SUBDOMAIN_NAMES = { "CYTOSOL", "ES", "MEMBRANE"};

//! \brief enum for subdomains (when using dune-multidomaingrid)
/* Important: Domains 0,..,nSubdomains-1 must correspond to the subdomain indices returned
 * by subGV.grid().domain(); -1 is a special value for the complete domain; further custom
 * values can be added with decreasing values starting with -2. This is especially useful
 * for identifying special gridfunctions which live on exotic subsets of entities (like my
 * beloved membrane interface GFs)
 */
enum GridDomains { DOMAIN_ALL = -1, DOMAIN_ELEC = 0, DOMAIN_MEMB = 1, DOMAIN_MEMB_INTERFACE = -2};
static const std::vector<std::string> MD_SUBDOMAINS = { "ELEC", "MEMB"};


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


#endif // DUNE_AX1_CONSTANTS_HH

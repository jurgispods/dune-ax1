#ifndef DUNE_AX1_ELECTROLYTE_HH
#define DUNE_AX1_ELECTROLYTE_HH

#include <valarray>
#include <vector>

#include <dune/ax1/common/constants.hh>

// Ion
template<class T>
class Ion
{
  
  public:
    Ion (T valence_, std::string name_, T relCon_=1.0)
      : valence(valence_), relCon(relCon_), name(name_)
    {}
  
    T getValence () const
    {
      return valence;
    }

    T getRelCon () const
    {
      return relCon;
    }

    std::string getName() const
    {
      return name;
    }

  private:
    T valence;
    T relCon; // relative concentration for stationary case
    T diffConst;
    std::string name;
};

// Solvent
template<class T>
class Solvent
{
private:
  T permittivity;
  
public:
  Solvent (T permittivity_)
    : permittivity(permittivity_)
  {}
  
  T getPermittivity () const { return permittivity; }
};

// Electrolyte
template<class T>
class Electrolyte
{
  public:
    Electrolyte (const T permittivity_, const T temperature_, const T stdCon_, const T lengthScale)
      : permittivity(permittivity_), temperature(temperature_), stdCon(stdCon_)
    {
      debyeLength = std::sqrt( 0.5 * con_eps0 * con_k * temperature / ( con_e * con_e * stdCon ) );
      //lengthConstantSqr = con_eps0 * con_k * temperature / ( con_e * con_e * stdCon );
      poissonConstant = con_e * con_e * stdCon * lengthScale * lengthScale / ( con_eps0 * con_k * temperature );
    }

    T getDebyeLength () const
    {
      return debyeLength;
    }

    T getPoissonConstant () const
    {
      return poissonConstant;
    }

    T getPermittivity () const
    {
      return permittivity;
    }

    void setPermittivity(T perm)
    {
      permittivity = perm;
    }

    T getTemperature () const
    {
      return temperature;
    }

    T getStdCon () const
    {
      return stdCon;
    }

    // add ion to electrolyte
    void addIon (Ion<T> ion)
    {
      ions.push_back(ion);
      con_diffWater.resize(ions.size());
    }

    // number of ion species
    int numOfSpecies () const
    {
      return ions.size();
    }

    // right hand side for the Poisson Boltzmann equation
    T rhsPoissonBoltzmann (const T phi) const
    {
      T sum = 0.0;
      for (int i=0; i<ions.size(); ++i)
      {
        sum = sum + ions[i].getValence() * ions[i].getRelCon() * exp(-ions[i].getValence() * phi);
      }
      return - 0.5 * sum / ( debyeLength * debyeLength );
    }

    // concentration of ion species for stationary case
    T getConcentration (const int& i, const T& phi) const
    {
      return stdCon * ions[i].getRelCon() * exp(-ions[i].getValence() * phi);
    }

    // get diffusion constant
    T getDiffConst ( const unsigned int ionSpecies ) const
    {
      assert(ionSpecies <= con_diffWater.size());
      return con_diffWater[ionSpecies];
    }

    void setDiffConst ( const unsigned int ionSpecies, T diffCoeff )
    {
      assert(ionSpecies <= con_diffWater.size());
      con_diffWater[ionSpecies] = diffCoeff;
    }
  
    // valence of ion species
    T getValence ( const unsigned int ionSpecies ) const
    {
      return ions[ionSpecies].getValence();
    }

    // name of ion species
    std::string getIonName ( const unsigned int ionSpecies ) const
    {
      return ions[ionSpecies].getName();
    }

    // charge density
    void addToChargeDensity(std::valarray<T>& chargeDensity,
                            const std::valarray<T>& concentrations,
                            const unsigned int ionSpecies)
    {
      chargeDensity += ions[ionSpecies].getValence() * concentrations;
    }

  private:
    T permittivity;
    std::vector<Ion<T> > ions;            // collection of ion species
    std::vector<T>       con_diffWater;   // corresponding diff coeffs for ions
    T temperature;
    T stdCon;                 // scaling concentration
    T debyeLength;
    T poissonConstant;
};

#endif // DUNE_AX1_ELECTROLYTE_HH

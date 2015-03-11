/*
 * power_parameters.hh
 *
 *  Created on: Aug 22, 2014
 *      Author: jpods
 */

#ifndef DUNE_AX1_POWER_PARAMETERS_HH
#define DUNE_AX1_POWER_PARAMETERS_HH

// TODO Derive from PowerGridFunction(Traits)? -> setTime() etc.
// Helper class to expose convection-diffusion parameter traits in array of parameters
template<class T, int k>
class PowerParameters :
  public std::vector<Dune::shared_ptr<T> >
{
  public:
    typedef std::vector<Dune::shared_ptr<T> > BaseT;
    typedef T ParamT;
    typedef typename T::Traits Traits;

    PowerParameters() : std::vector<Dune::shared_ptr<T> >()
    {}

    void setTime(double time)
    {
      for(int i=0; i<this->size(); ++i)
      {
        this->operator[](i)->setTime(time);
      }
    }

    void preAssembly(bool implicit)
    {
//      for(int i=0; i<this->size(); ++i)
//      {
//        this->operator[](i)->preAssembly();
//      }
      // Only one parameter class need to update the (shared) membrane flux GFs!
      this->operator[](0)->preAssembly(implicit);
    }

};

#endif /* DUNE_AX1_POWER_PARAMETERS_HH */

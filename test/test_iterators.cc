#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <stdlib.h>
#include <typeinfo>

#include <dune/common/typetraits.hh>

#include <dune/ax1/common/constants.hh>
//#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_interfacevector.hh>

struct Blob
{
  Blob(int d_)
  : d(d_)
  {}

  int f() const
  {
    return 3;
  }

  int d;
};


//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  typedef Blob T;
  std::vector<T> v1;
  for(int i=0; i<10; i++)
  {
    T t(i);
    v1.push_back(t);
  }

  typedef std::vector<T*> Container;
  Container pVec;
  for(int i=0; i<5; i++)
  {
    T* pInt = &v1[v1.size()-1-i];
    pVec.push_back(pInt);
    //printf("%d\n", (*mybegin).f());
    std::cout << "pVec[" << i << "] = " << (*pVec[i]).d << std::endl;
    std::cout << "pVec[" << i << "] = " << pVec[i]->d << std::endl;
  }

  typedef Ax1InterfaceIterator<Container,T*,T> MyIterator;
  MyIterator mybegin(pVec, 0);
  printf("%d\n", (*mybegin).f());
  //std::cout << "mybegin = " << (*mybegin) << std::endl;


  MyIterator myend(pVec, pVec.size());

  std::cout << "mybegin == myend ? " << (mybegin == myend) << std::endl;

  for(MyIterator mit = mybegin; mit != myend; mit++)
  {
    std::cout << "mit = " << (*mit).d << std::endl;
    printf("%d\n", mit->f());
  }

  const std::vector<T> v1_const(v1);



  typedef const std::vector<T*> ConstContainer;
  ConstContainer pVec_const(pVec);
//  for(int i=0; i<5; i++)
//  {
//    const int* pInt = &v1_const[v1_const.size()-1-i];
//    pVec_const.push_back(pInt);
//  }
//
//  typedef Dune::GenericIterator<ConstContainer,const int*,const int*&> DuneConstIterator;
//  DuneConstIterator dunebegin_const(pVec_const, 0);
//  DuneConstIterator duneend_const(pVec_const, pVec_const.size());
// for(DuneConstIterator dit = dunebegin_const; dit != duneend_const; dit++)
// {
//   std::cout << "dit = " << (*dit) << std::endl;
// }


  typedef Ax1InterfaceIterator<ConstContainer,T*,const T> MyConstIterator;
  MyConstIterator mybegin_const(pVec_const, 0);
  MyConstIterator myend_const(pVec_const, pVec_const.size());
  for(MyConstIterator mit = mybegin_const; mit != myend_const; ++mit)
  {
    std::cout << "mit = " << mit->d << std::endl;
  }


  return 0;

}




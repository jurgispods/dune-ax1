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

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/typetraits.hh>

#include <dune/ax1/common/tools.hh>
#include <dune/ax1/common/ax1_c++11_geekstuff.hh>

#define HAS_MEM_FUNC(func, name)                                        \
    template<typename T, typename Sign>                                 \
    struct name {                                                       \
        typedef char yes[1];                                            \
        typedef char no [2];                                            \
        template <typename U, U> struct type_check;                     \
        template <typename C> static yes &chk(type_check<Sign, &C::func> *); \
        template <typename C> static no  &chk(...);                    \
        static bool const value = sizeof(chk<T>(0)) == sizeof(yes);     \
    }

#define HAS_TYPEDEF(typed, name)                                        \
    template<typename T>                                 \
    struct name {                                                       \
        typedef char yes[1];                                            \
        typedef char no [2];                                            \
        template <typename U> struct typedef_check;                     \
        template <typename C> static yes &chk(typedef_check<typename C::typed> *); \
        template <typename C> static no  &chk(...);                    \
        static bool const value = sizeof(chk<T>(0)) == sizeof(yes);     \
    }

namespace detail
{
  template<template<typename,typename> class V>
  struct MyStruct
  {
     typedef int type;
  };

   template<>
   struct MyStruct<std::vector>
   {
     typedef double type;
   };
}

template<class T>
class MyClass
{
  public:
    typedef int type;

    void info()
    {
      typename detail::MyStruct<std::vector>::type bla = 0;
      std::cout << "Type of bla: " << Tools::getTypeName(bla) << std::endl;
      std::cout << "bla = " << bla << std::endl;

      typename detail::MyStruct<std::list>::type bla2 = 0;
      std::cout << "Type of bla2: " << Tools::getTypeName(bla2) << std::endl;
      std::cout << "bla2 = " << bla2 << std::endl;
    }

};


template<class T>
class MyDerivedClass : public MyClass<T>
{
public:

//
};


template <typename T, int (T::*) ()> struct enable { typedef T type; };
template <typename T> typename enable<T, &T::i>::type bla (T&);

struct A
{
    typedef int WurstBrot;

    void i()
    {
      std::cout << "hallo" << std::endl;
    }
};

struct B
{
    int i()
    {
      std::cout << "hallo" << std::endl;
      return 0;
    }

    std::string toString()
    {
      return std::string("wurst");
    }
};

template<typename T>
struct C
{
  typedef int WurstBrot;
  typedef std::vector<T> Vector;

  T t;

  void hurz(T t)
  {
    std::cout << "hurz!" << std::endl;
  }

  void lurz()
  {
    std::cout << "lurz!" << std::endl;
  }
};

// Typedef inheritance does not work when base class is a
// template!
template<typename T>
struct D: public C<T>
{


  //WurstBrot w; // does not work
  //Vector v; // does not work
  void hallo()
  {
    this->t = 0;
    this->hurz(this->t);
    this->lurz();
  }
};

struct E: public C<int>
{
  WurstBrot w;
};


HAS_MEM_FUNC(toString, has_to_string);

template<typename T> void
doSomething() {

  typename std::string(T::*hurz)();
  std::cout << Tools::getTypeName(hurz) << std::endl;

   if(has_to_string<T, std::string(T::*)()>::value) {

      std::cout << " has function toString()!";
   } else {
     std::cout << " does not have function toString()!";
   }
   std::cout << std::endl;
}

// TODO Das geht natuerlich nicht, ein typedef kann
// nicht als Fumntions-Pointer gefunden werden. Das Marko
// fuer die Memberfunktionen kann also nicht genutzt werden,
// also alternatives Makro schreiben!
HAS_TYPEDEF(WurstBrot, has_wurstbrot);
template<typename T> void
doSomethingDifferent() {

   if(has_wurstbrot<T>::value) {

      std::cout << " has typedef WurstBrot!";
   } else {
     std::cout << " does not have typedef WurstBrot!";
   }
   std::cout << std::endl;
}

template<typename> struct void_ { typedef void type; };

template<typename T, typename = void>  // Line 1
struct is_class { static bool const value = false; };

// This will get preference if it exists
// (This is case when 'int T::*', i.e. any member pointer, exists)
template<typename T>
struct is_class<T, typename void_<int T::*>::type> { // Line 2
  static bool const value = true;
};

template<typename T, typename = void>  // Line 1
struct has_typedef {
  static bool const value = false;
};

template<typename T>
struct has_typedef<T, typename void_<typename T::WurstBrot>::type> { // Line 2
  static bool const value = true;
};

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  /*

  typedef MyClass<int> T1;
  typedef MyDerivedClass<int> T2;
  typedef MyDerivedClass<float> T3;

  int hurz = Dune::IsBaseOf<T1,T2>::value;
  int hurz2 = Dune::IsBaseOf<T2,T1>::value;
  int hurz3 = Dune::IsBaseOf<T1,T3>::value;
  //int hurz4 = Dune::IsBaseOf<MyClass,MyDerivedClass>::value;

  std::cout << "hurz = " << hurz << std::endl;
  std::cout << "hurz2 = " << hurz2 << std::endl;
  std::cout << "hurz3 = " << hurz3 << std::endl;
  //std::cout << "hurz4 = " << hurz4 << std::endl;

  MyClass<int> asd1;
  asd1.info();
  MyDerivedClass<int> asd2;
  //asd2.hurz();


  A a; a.i();
  B b; b.i();
  //bla(a);
  //bla(b);
  C<int> c;
  D<int> d;
  d.hallo();
  //std::cout << Tools::getTypeName(d.v) << std::endl;
  E e;
  std::cout << Tools::getTypeName(e.w) << std::endl;


  doSomething<A>();
  doSomething<B>();

  doSomethingDifferent<A>();
  doSomethingDifferent<B>();

  std::cout << is_class<int>::value << std::endl;
  std::cout << is_class<A>::value << std::endl;

  std::cout << has_typedef<A>::value << std::endl;
  std::cout << has_typedef<B>::value << std::endl;


//  MyClass<float> asd2;
//  asd2.info();
//  MyClass<double> asd3;
//  asd3.info();



  std::string kacke = "Kacke.";
  std::cout << kacke << std::endl;

  */

  auto a = std::make_tuple(5, "Hello", -0.1);
  std::cout << a << std::endl; // prints: (5, "Hello", -0.1)

  return 0;

}




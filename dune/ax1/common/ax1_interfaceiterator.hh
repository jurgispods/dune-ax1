/*
 * ax1_interfacevector.hh
 *
 *  Created on: Jun 6, 2013
 *      Author: jpods
 */

#ifndef DUNE_AX1_INTERFACEVECTOR_HH
#define DUNE_AX1_INTERFACEVECTOR_HH


#include <dune/common/iteratorfacades.hh>
#include<dune/common/typetraits.hh>
#include <cassert>

/**
 * @brief Devil in disguise!
 *
 * Custom iterator specifically designed for underlying containers of Iterators (this one was written for
 * Dune::IntersectionIterators). Dereferencing the iterator returns a reference of type R, which is obtained
 * by dereferencing the underlying RIterator, in the
 * specific case of an IntersectionIterator this would correspond to an Intersection.
 *
 * \tparam C The underlying container with values of type RIterator
 * \tparam RIterator The element type of the container, must be pointer-like, i.e. support a dereferencing operator!
 * \tparam R The return type of the dereferencing operator of the RIterator class. This is the value that is returned
 * when calling the dereferencing operator!
 *
 * Please note: Due to my lack of TMP knowledge, this class is not written in a very generic way. Here is an example
 * of know to define mutable and const iterators for the case of RIterator = IntersectionIterator, R = Intersection:
 *
 * \code
 *
 * typedef std::vector<IntersectionIterator> IVector;
 * //Note: Just put in 'pure' type for template param R, reference qualifier is added in class!
 * typedef Ax1InterfaceIterator<IVector, IntersectionIterator, Intersection> MyIterator;
 * //Note: It is enough to make C and R const for a const iterator, RIterator has to be mutable!
 * typedef Ax1InterfaceIterator<const IVector, IntersectionIterator, const Intersection> MyConstIterator;
 *
 * IVector v = setupVector();
 * MyIterator mybegin(v,0);
 * MyIterator mend(v,v.size());
 *
 * for(MyIterator mit = mybegin; mit != mend; ++mit)
 * {
 *   std::cout << mit->geometry().center() << std::endl;
 * }
 *
 * \endcode
 *
 */
template<class C, class RIterator, class R>
class Ax1InterfaceIterator :
    public std::iterator< std::forward_iterator_tag,
              RIterator, // std::iterator needs mutable value type
              std::ptrdiff_t,
              R*,
              R&>
{

  typedef Ax1InterfaceIterator<C,RIterator,R> This;

public:

  typedef C Container;
  typedef RIterator Value;
  typedef std::ptrdiff_t DifferenceType;
  typedef R& Reference;
  typedef R* Pointer;

  // Constructors needed by the base iterators
  Ax1InterfaceIterator()
  : container_(0),
    position_(0),
    piType(Dune::PartitionIteratorType::All_Partition)
  {
    init();
  }

  /**
   * @brief Constructor
   * @param cont Reference to the container we are an iterator for
   * @param pos The postion the iterator will be positioned to
   * (e.g. 0 for an iterator returned by Container::begin() or
   * the sizeof the container for an iterator returned by Container::end()
   */
  Ax1InterfaceIterator(Container& cont, DifferenceType pos,
      Dune::PartitionIteratorType piType_ = Dune::PartitionIteratorType::All_Partition)
    : container_(&cont),
      position_(pos),
      piType(piType_)
  {
    init();
  }

  /**
   * @brief Copy constructor
   *
   * This is somehow hard to understand, therefore play with the cases:
   * 1. if we are mutable this is the only valid copy constructor, as the argument is a mutable iterator
   * 2. if we are a const iterator the argument is a mutable iterator => This is the needed conversion to initialize a const iterator from a mutable one.
   */
  Ax1InterfaceIterator(const This& other):
    container_(other.container_),
    position_(other.position_),
    piType(other.piType)
  {
    init();
  }


  // Methods needed by the forward iterator
  bool equals(const This & other) const
  {
    return position_ == other.position_ && container_ == other.container_;
  }

  Reference dereference() const{
    return *container_->operator[](position_);
  }

  Pointer operator->() const
  {
    return &(this->dereference());
  }

  /** @brief Dereferencing operator. */
  Reference operator*() const
  {
    return dereference();
  }

  bool operator==(const This& rhs)
  {
    return this->equals(rhs);
  }

  bool operator!=(const This& rhs)
    {
      return !this->equals(rhs);
    }

  void increment(){
    ++position_;
  }

  /** @brief Preincrement operator. */
  This& operator++()
  {
    this->increment();
    // Increment until the desired partition type is found!
    while(position_ < container_->size() && !hasMatchingPartitionType())
    {
      //debug_jochen << "ptype " << dereference().inside()->partitionType() << " of is @" << dereference().geometry().center()
      //    << " does not match, incrementing!" << std::endl;
      this->increment();
    }
    return *this;
  }

  /** @brief Postincrement operator. */
  This operator++(int)
  {
    This tmp(*this);
    this->operator++();
    return tmp;
  }


private:

  void init()
  {
    // Move position to the first occurrence of a matching partition type
    if(position_ < container_->size() && !hasMatchingPartitionType())
    {
      this->operator++();
    }
  }

  bool hasMatchingPartitionType()
  {
    switch(piType)
    {
      case Dune::PartitionIteratorType::All_Partition:
        return true;
        break;
      case Dune::PartitionIteratorType::Interior_Partition:
        if(dereference().inside()->partitionType() == Dune::PartitionType::InteriorEntity)
          return true;
        else
          return false;
       break;
      default:
        DUNE_THROW(Dune::NotImplemented, "Partition iterator type '" << piType
            << "' not implemented yet, only All_Partition and Interior_Partition are available!");

    }
  }

  Container *container_;
  DifferenceType position_;
  Dune::PartitionIteratorType piType;
};





//template<class T,class R>
//class InterfaceVector : public std::vector<T>
//{
//public:
//  typedef Ax1InterfaceIterator<InterfaceVector<T,R>,T,R> iterator;
//  typedef Ax1InterfaceIterator<const InterfaceVector<T,R>,const T,const R&> const_iterator;
//
//
//  iterator begin(){
//    return iterator(*this, 0);
//  }
//
//  const_iterator begin() const{
//    return const_iterator(*this, 0);
//  }
//
//  iterator end(){
//    return iterator(*this, size());
//  }
//
//  const_iterator end() const{
//    return const_iterator(*this, size());
//  }
//
//};

#endif /* DUNE_AX1_INTERFACEVECTOR_HH */

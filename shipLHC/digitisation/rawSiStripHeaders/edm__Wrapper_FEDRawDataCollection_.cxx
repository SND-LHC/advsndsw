#include "edm__Wrapper_FEDRawDataCollection_.h"
using namespace std;

edm::Wrapper<FEDRawDataCollection>::Wrapper() {
}
edm::Wrapper<FEDRawDataCollection> &edm::Wrapper<FEDRawDataCollection>::operator=(const Wrapper & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   edm::WrapperBase::operator=(const_cast<Wrapper &>( rhs ));
   present = (const_cast<Wrapper &>( rhs ).present);
   obj = (const_cast<Wrapper &>( rhs ).obj);
   return *this;
}
edm::Wrapper<FEDRawDataCollection>::Wrapper(const Wrapper & rhs)
   : edm::WrapperBase(const_cast<Wrapper &>( rhs ))
   , obj(const_cast<Wrapper &>( rhs ).obj)
   , present(const_cast<Wrapper &>( rhs ).present)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
edm::Wrapper<FEDRawDataCollection>::~Wrapper() {
}
#endif // edm__Wrapper_FEDRawDataCollection__cxx
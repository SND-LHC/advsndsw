#include "edm__WrapperBase.h"
using namespace std;

edm::WrapperBase::WrapperBase() {
}
edm::WrapperBase &edm::WrapperBase::operator=(const WrapperBase & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   edm::ViewTypeChecker::operator=(const_cast<WrapperBase &>( rhs ));
   return *this;
}
edm::WrapperBase::WrapperBase(const WrapperBase & rhs)
   : edm::ViewTypeChecker(const_cast<WrapperBase &>( rhs ))
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
edm::WrapperBase::~WrapperBase() {
}
#endif // edm__WrapperBase_cxx

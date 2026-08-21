#include "edm__ViewTypeChecker.h"
using namespace std;

edm::ViewTypeChecker::ViewTypeChecker() {
}
edm::ViewTypeChecker &edm::ViewTypeChecker::operator=(const ViewTypeChecker & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   return *this;
}
edm::ViewTypeChecker::ViewTypeChecker(const ViewTypeChecker & rhs)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
edm::ViewTypeChecker::~ViewTypeChecker() {
}
#endif // edm__ViewTypeChecker_cxx
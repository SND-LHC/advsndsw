#include "edm__DoNotRecordParents.h"
using namespace std;

edm::DoNotRecordParents::DoNotRecordParents() {
}
edm::DoNotRecordParents &edm::DoNotRecordParents::operator=(const DoNotRecordParents & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   return *this;
}
edm::DoNotRecordParents::DoNotRecordParents(const DoNotRecordParents & rhs)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
}
edm::DoNotRecordParents::~DoNotRecordParents() {
}
#endif // edm__DoNotRecordParents_cxx

#include "FEDRawDataCollection.h"
using namespace std;

FEDRawDataCollection::FEDRawDataCollection() {
}
FEDRawDataCollection &FEDRawDataCollection::operator=(const FEDRawDataCollection & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   edm::DoNotRecordParents::operator=(const_cast<FEDRawDataCollection &>( rhs ));
   data_ = (const_cast<FEDRawDataCollection &>( rhs ).data_);
   FEDRawDataCollection &modrhs = const_cast<FEDRawDataCollection &>( rhs );
   modrhs.data_.clear();
   return *this;
}
FEDRawDataCollection::FEDRawDataCollection(const FEDRawDataCollection & rhs)
   : edm::DoNotRecordParents(const_cast<FEDRawDataCollection &>( rhs ))
   , data_(const_cast<FEDRawDataCollection &>( rhs ).data_)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   FEDRawDataCollection &modrhs = const_cast<FEDRawDataCollection &>( rhs );
   modrhs.data_.clear();
}
FEDRawDataCollection::~FEDRawDataCollection() {
}
#endif // FEDRawDataCollection_cxx
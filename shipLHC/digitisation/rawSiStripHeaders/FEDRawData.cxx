#include "FEDRawData.h"
using namespace std;

FEDRawData::FEDRawData() {
}
FEDRawData &FEDRawData::operator=(const FEDRawData & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   data_ = (const_cast<FEDRawData &>( rhs ).data_);
   FEDRawData &modrhs = const_cast<FEDRawData &>( rhs );
   modrhs.data_.clear();
   return *this;
}
FEDRawData::FEDRawData(const FEDRawData & rhs)
   : data_(const_cast<FEDRawData &>( rhs ).data_)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   FEDRawData &modrhs = const_cast<FEDRawData &>( rhs );
   modrhs.data_.clear();
}
FEDRawData::~FEDRawData() {
}
#endif // FEDRawData_cxx
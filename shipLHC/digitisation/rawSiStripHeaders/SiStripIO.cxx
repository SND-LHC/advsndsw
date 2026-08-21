using namespace std;
#include "SiStripIO.h"

#ifndef edm__Wrapper_FEDRawDataCollection__cxx
#define edm__Wrapper_FEDRawDataCollection__cxx
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

#ifndef edm__WrapperBase_cxx
#define edm__WrapperBase_cxx
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

#ifndef edm__ViewTypeChecker_cxx
#define edm__ViewTypeChecker_cxx
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

#ifndef FEDRawDataCollection_cxx
#define FEDRawDataCollection_cxx
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

#ifndef edm__DoNotRecordParents_cxx
#define edm__DoNotRecordParents_cxx
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

#ifndef FEDRawData_cxx
#define FEDRawData_cxx
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
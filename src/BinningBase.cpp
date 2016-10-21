#include "BinningBase.h"

BinningBase::BinningBase() :
  _dimension  (0 ),
  _axisNames  (0 ),
  _binningType("")   
{


}

bool BinningBase::isSameBinningType( const BinningBase& other ) const{
  return (other._binningType == _binningType);
}

void BinningBase::setBinningType(TString binningType){
  _binningType = binningType;
}


void      BinningBase::setNames( HyperName names ){ 
  _axisNames = names; 
}  

HyperName BinningBase::getNames() const{
  return _axisNames;
}       


const int& BinningBase::getDimension () const{ 
  return _dimension; 
}  

void BinningBase::setDimension (int dimension){
  if (_dimension == 0){
    _dimension       = dimension;
    _axisNames       = HyperName  (dimension);
  }
}


double BinningBase::getMin(int dimension) const{
  return getLimits().getLowCorner().at(dimension);
}

double BinningBase::getMax(int dimension) const{
  return getLimits().getHighCorner().at(dimension);
}
  
TString BinningBase::getBinningType() const{
  return _binningType;
}

bool BinningBase::isDiskResident() const{
  return false;
}
TString BinningBase::filename() const{
  return "";
}

BinningBase::~BinningBase(){

}
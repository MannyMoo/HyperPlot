#include "HyperBinningMakers.h"


/***************************************************************************
                         HyperBinningMakerSmart
***************************************************************************/


HyperBinningMakerSmart::HyperBinningMakerSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerSmart() Constructor";  
  //makeBinning(startingDim);
}

void HyperBinningMakerSmart::makeBinning(){
   
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim; 
  if (splitDim >= dimension) splitDim = 0; 
  
  int nBins = 0;
  int unchanged = 0;  

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim) << std::endl;

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0; 
    if (unchanged > dimension) break;
    nBins = getNumBins(); 
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<< std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<< std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Smart binning algorithm complete"<< std::endl;

}

HyperBinningMakerSmart::~HyperBinningMakerSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerSmart() Constructor";  
}

/***************************************************************************
                         HyperBinningMakerMint
***************************************************************************/


HyperBinningMakerMint::HyperBinningMakerMint(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerMint() Constructor"<<std::endl; 
  //makeBinning(startingDim);
}

void HyperBinningMakerMint::makeBinning(){
  
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim;
  if (splitDim >= dimension) splitDim = 0;

  int nBins = 0;
  int unchanged = 0;

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim) <<std::endl;

  while ( splitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins" << std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim) <<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Mint binning algorithm complete " <<std::endl;

}
 
HyperBinningMakerMint::~HyperBinningMakerMint(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerMint() Constructor" <<std::endl; 
}


/***************************************************************************
                         HyperBinningMakerMultiSmart
***************************************************************************/


HyperBinningMakerMultiSmart::HyperBinningMakerMultiSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerMultiSmart() Constructor";  
  //makeBinning(startingDim);
}

void HyperBinningMakerMultiSmart::makeBinning(){
   
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim; 
  if (splitDim >= dimension) splitDim = 0; 
  
  int nBins = 0;
  int unchanged = 0;  

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim) << std::endl;

  while (smartMultiSplitAll(_binningDimensions.at(splitDim)) != -1){
    finishedIteration();
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0; 
    if (unchanged > dimension) break;
    nBins = getNumBins(); 
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<< std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<< std::endl;
  }
  
  INFO_LOG << "Gone as far as possible with MultiSplit - now trying SmartSplit " << std::endl;
  unchanged = 0;  

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != -1){
    finishedIteration();
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0; 
    if (unchanged > dimension) break;
    nBins = getNumBins(); 
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<< std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<< std::endl;
  }


  if (s_printBinning == true) INFO_LOG << "Smart binning algorithm complete"<< std::endl;

}

HyperBinningMakerMultiSmart::~HyperBinningMakerMultiSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerMultiSmart() Constructor";  
}


/***************************************************************************
                         HyperBinningMakerMintSmart
***************************************************************************/


HyperBinningMakerMintSmart::HyperBinningMakerMintSmart(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerMintSmart() Constructor"<<std::endl; 
  //makeBinning(startingDim);
}

void HyperBinningMakerMintSmart::makeBinning(){
  
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim;
  if (splitDim >= dimension) splitDim = 0;

  int nBins = 0;
  int unchanged = 0;

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;

  while ( splitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;
  }

  nBins = 0;
  unchanged = 0;

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged > dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Mint binning algorithm complete "<<std::endl;

}

HyperBinningMakerMintSmart::~HyperBinningMakerMintSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerMintSmart() Constructor"<<std::endl;  
}


/***************************************************************************
                         HyperBinningMakerSmartRandomise
***************************************************************************/


HyperBinningMakerSmartRandomise::HyperBinningMakerSmartRandomise(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerRandomise() Constructor"<<std::endl;  
}

void HyperBinningMakerSmartRandomise::makeBinning(){
  
  if (s_printBinning == true) INFO_LOG << "Splitting all bins in random dimensions"<<std::endl;
  int dimension = _binningDimensions.size();  
  
  int nBins = 0;
  int unchanged = 0;

  while (smartSplitAllRandomise() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Random binning algorithm complete "<<std::endl;


}

HyperBinningMakerSmartRandomise::~HyperBinningMakerSmartRandomise(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerRandomise() Constructor";  
}


/***************************************************************************
                         HyperBinningMakerMintRandomise
***************************************************************************/


HyperBinningMakerMintRandomise::HyperBinningMakerMintRandomise(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerRandomise() Constructor"<<std::endl;  
}

void HyperBinningMakerMintRandomise::makeBinning(){
  
  if (s_printBinning == true) INFO_LOG << "Splitting all bins in random dimensions"<<std::endl;
  int dimension = _binningDimensions.size();  

  int nBins = 0;
  int unchanged = 0;

  while (splitAllRandomise() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Random binning algorithm complete "<<std::endl;


}

HyperBinningMakerMintRandomise::~HyperBinningMakerMintRandomise(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerRandomise() Constructor"<<std::endl;  
}



/***************************************************************************
                         HyperBinningMakerLikelihood
***************************************************************************/


HyperBinningMakerLikelihood::HyperBinningMakerLikelihood(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerLikelihood() Constructor"<<std::endl; 
}

void HyperBinningMakerLikelihood::makeBinning(){
  
  int dimension = _binningDimensions.size();  
  
  int nBins = 0;
  int unchanged = 0;

  while (likelihoodSplitAll() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "likelihood binning algorithm complete "<<std::endl;

}

HyperBinningMakerLikelihood::~HyperBinningMakerLikelihood(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerLikelihood() Constructor"<<std::endl; 
}


/***************************************************************************
                         HyperBinningMakerSmartLikelihood
***************************************************************************/


HyperBinningMakerSmartLikelihood::HyperBinningMakerSmartLikelihood(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerSmartLikelihood() Constructor"<<std::endl;  
}

void HyperBinningMakerSmartLikelihood::makeBinning(){
  
  int dimension = _binningDimensions.size();  
  
  int nBins = 0;
  int unchanged = 0;

  while (smartLikelihoodSplitAll() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins" << std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "likelihood binning algorithm complete " << std::endl;

}

HyperBinningMakerSmartLikelihood::~HyperBinningMakerSmartLikelihood(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerSmartLikelihood() Constructor"<<std::endl; 
}





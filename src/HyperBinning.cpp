#include "HyperBinning.h"


///The only constructor
HyperBinning::HyperBinning() :
  _changed(true),
  _averageBinWidth(getDimension()),
  _minmax( HyperPoint(getDimension()), HyperPoint(getDimension()) )
{
  setBinningType("HyperBinning");
  WELCOME_LOG << "Hello from the HyperBinning() Constructor";
}



///Set the dimension of the HyperBinning. This can only be 
///called once, when it is known what dimesnion it is.
void HyperBinning::setDimension(int dim){
  
  if (getDimension() == 0){
    BinningBase::setDimension(dim);
    _averageBinWidth = HyperPoint ( getDimension(), 1.0);
    _minmax          = HyperCuboid( getDimension(), 0.0, 1.0 );
  }

}


/** 
    This is used to get the bin number that the HyperPoint falls into.

    This is done by looping over the HyperVolumes in order until the HyperPoint
    is contained within one of them. When one of these is found, check to see if there
    are any 'linked bins' - if not, this is a true bin, so job done. If there are
    linked bins, proceed to see which of the linked bins the HyperPoint falls into. 
    Continue this process until you reach a HyperVolume with no linked bins.

    The BinNumber is found from the HyperVolume number using the lookup vector _binNum.
             
    ~~~ {.cpp}
                  0   1   2   3   4   5   6   7   8
     _binNum = { -1, -1, -1, -1,  2,  3,  4,  0   1 }


           HyperVolume Numbers 
    
     |-------------0-------------| 
    
     |------1------|------2------| 
    
     |--3---|---4--|---5---|--6--|
    
     |-7-|-8| 
    
               Bin Numbers
    
     | 0 | 1|   2  |   3   |  4  |
    
    ~~~
*/
int HyperBinning::getBinNum(const HyperPoint& coords) const{
  
  if (_changed == true) updateCash();

  //First check if the HyperPoint is in the HyperCuboid _minmax that
  //surrounds all the bins.

  if (_minmax.inVolume(coords) == 0) return -1;
  
  int nPrimVols = getNumPrimaryVolumes();
  
  if ( nPrimVols == 0){

    //loop over all the bins until one is found that contains the event. If a bin
    // is found that contains the event, check if it has any linked bins. Linked bins
    // aren't real bins... just used to speed up sorting later.

    int volumeNumber = -1;
  
    for (int i = 0; i < getNumHyperVolumes(); i++){
      bool inVol = getHyperVolume(i).inVolume(coords);
      if (inVol == 1) { volumeNumber = i; break; }
    }
     
    if (volumeNumber == -1) return -1;
  
    if ( getLinkedHyperVolumes(volumeNumber).size() > 0 ) volumeNumber = followBinLinks(coords, volumeNumber);
  
    return _binNum.at(volumeNumber);
  }


  int primaryVolumeNumber = -1;

  for (int i = 0; i < nPrimVols; i++){
    int thisVolNum = getPrimaryVolumeNumber(i);
    bool inVol = getHyperVolume(thisVolNum).inVolume(coords);
    if (inVol == 1) { primaryVolumeNumber = thisVolNum; break; }
  }
  
  int volumeNumber = -1;

  if ( getLinkedHyperVolumes(primaryVolumeNumber).size() > 0 ) {
    volumeNumber = followBinLinks(coords, primaryVolumeNumber);
  }
  else{
    ERROR_LOG << "This primary volume has NO links. Not what I expect!!" << std::endl;
    return _binNum.at(primaryVolumeNumber);
  }
  
  return _binNum.at(volumeNumber);

}



bool HyperBinning::isPrimaryVolume(int volumeNumber) const{
  
  int nPrimVols = getNumPrimaryVolumes();

  for (int i = 0; i < nPrimVols; i++){
    if (getPrimaryVolumeNumber(i) == volumeNumber) return true;
  }

  return false;

}





///Used to follow the bin hierarchy. Give it a HyperPoint, and the number of a 
///HyperVolume (that has links) that the HyperPoint falls into. 
///
///
int HyperBinning::followBinLinks(const HyperPoint& coords, int motherVolumeNumber) const{
  
  //find the linked volumes
  const std::vector<int> linkedVolumes = getLinkedHyperVolumes(motherVolumeNumber);
  
  int volumeNumber = -1;
  
  //see if the coords falls into any of the linked volumes (it should if there are no bugs)
  for (unsigned i = 0; i < linkedVolumes.size(); i++){
    int daughBinNum = linkedVolumes.at(i);
    bool inVol = getHyperVolume(daughBinNum).inVolume(coords);
    if (inVol == 1) { volumeNumber = daughBinNum; break; }
  }
  
  if (volumeNumber == -1) {
    ERROR_LOG << "The trail of linked bins has gone cold!";
    return -1;
  }
  
  //now have volumeNumber which contains the next bin in the hierarchy.
  // if this is linked to more bins, keep following the trail!
  if ( getLinkedHyperVolumes(volumeNumber).size() > 0 ) volumeNumber = followBinLinks(coords, volumeNumber);
  
  //if not, we have made it to the end. Return the volume number!
  return volumeNumber;

}



///Get number of bins (this is NOT the number of
///HyperVolumes!!! - see the class description for more details)
int HyperBinning::getNumBins() const{
  
  //std::cout << "getNumBins" << std::endl;
  if (_changed == true) updateCash();
  
  return _hyperVolumeNumFromBinNum.size();

}

///Get the HyperVolume assosiated with bin number i. This just uses 
///the _hyperVolumeNumFromBinNum variable to find the HyperVolume number
///from the bin number, then returns that HyperVolume.
HyperVolume HyperBinning::getBinHyperVolume(int binNumber) const{
  
  //std::cout << "getBinHyperVolume" << std::endl;
  if (_changed == true) updateCash();
  return getHyperVolume( _hyperVolumeNumFromBinNum.at(binNumber) );

}

///Get the bin number assosiated with a given HyperVolume number. 
///If this returns -1, it means that the HyperVolume in question
///is not a bin, but part of the binning hierarchy.
int HyperBinning::getBinNum(int volumeNumber) const{
  if (_changed == true) updateCash();
  return _binNum.at(volumeNumber);
}

/// get the HyperVolume Number from the bin number
///
int HyperBinning::getHyperVolumeNumber(int binNumber) const{
  if (_changed == true) updateCash(); 
  return _hyperVolumeNumFromBinNum.at(binNumber);
}

///Update the cash which includes the  mutable member variables
///_binNum, _hyperVolumeNumFromBinNum, _averageBinWidth,
/// and _minmax.
void HyperBinning::updateCash() const{
  
  //INFO_LOG << "Updating cache" << std::endl;

  //note the ordering of theses is important...
  //possilbe infinite loops

  _changed = false;

  updateBinNumbering();
  updateAverageBinWidth();
  updateMinMax();

}

///Update the member variables _binNum and _hyperVolumeNumFromBinNum.
///Will usually be called from updateCash()
void HyperBinning::updateBinNumbering() const{
  
  //first fill all the bin numbers with -1
  int nVolumes = getNumHyperVolumes();
  _binNum = std::vector<int>(nVolumes, -1);
  

  //if a HyperVolume has any linked HyperVolumes,
  //then set its bin number to count.
  int count = 0;
  for (int i = 0; i < getNumHyperVolumes(); i++){
    if ( getLinkedHyperVolumes(i).size() == 0 ) {
      _binNum.at(i) = count;
      count++;
    }
  }  

  //now we know how many bins there are, make the 
  // _hyperVolumeNumFromBinNum vector
  int nBins = count;
  _hyperVolumeNumFromBinNum = std::vector<int>(nBins,-1);

  //fill the vector
  for (int i = 0; i < getNumHyperVolumes(); i++){
    if ( _binNum.at(i) != -1 ) {
      _hyperVolumeNumFromBinNum.at( _binNum.at(i) ) = i;
    }
  }    

 
}

///return the limits of the binning.
///This value is cashed for speed - when the binning changes the cashe will
///automatically be updated.
HyperCuboid HyperBinning::getLimits() const{
  if (_changed == true) updateCash();   
  return _minmax;
}



///update the _averageBinWidth HyperPoint. 
///Will usually be called from updateCash()
void HyperBinning::updateAverageBinWidth() const{
  
  int dim = getDimension();

  HyperPoint averageWidth(dim);

  for (int i = 0; i < getNumBins(); i++){
    for (int j = 0; j < dim; j++) {
      double min = getBinHyperVolume(i).getMin(j);
      double max = getBinHyperVolume(i).getMax(j);
      averageWidth.at(j) += (max - min);
    }  
  }    

  _averageBinWidth = averageWidth/(double)getNumBins();
  
  //_averageBinWidth.print();
}


///Update the miniumum and maximum values, _minmax, 
///in the cashe. Will usually be called from updateCash().
void HyperBinning::updateMinMax() const{
  
  int dim = getDimension();

  HyperPoint min(dim);
  HyperPoint max(dim);
  
  for (int d = 0; d < dim; d++){
    min.at(d) = getHyperVolume(0).getMin(d);
    max.at(d) = getHyperVolume(0).getMax(d);
  }

  for(int i = 1; i < getNumHyperVolumes(); i++){
    HyperVolume thisVol = getHyperVolume(i);
    for (int d = 0; d < dim; d++){
      if (min.at(d) > thisVol.getMin(d)) min.at(d) = thisVol.getMin(d);
      if (max.at(d) < thisVol.getMax(d)) max.at(d) = thisVol.getMax(d);
    }
  }

  _minmax = HyperCuboid(min, max);

}


///get the average bin width HyperPoint (average bin width in each dimension).
///This value is cashed for speed - when the binning changes the cashe will
///automatically be updated.
HyperPoint HyperBinning::getAverageBinWidth() const{
  //std::cout << "getAverageBinWidth" << std::endl;
  if (_changed == true) updateCash();
  return _averageBinWidth;  

}


void HyperBinning::reserveCapacity(int nElements){

  _binNum                  .reserve(nElements);
  _hyperVolumeNumFromBinNum.reserve(nElements);

}



/// This function will merge the two binnings. It assumes that the first
///
void HyperBinning::mergeBinnings( const BinningBase& other ){
  
  //INFO_LOG << "Starting HyperBinning::mergeBinnings" << std::endl;

  if ( other.getBinningType().Contains("HyperBinning") == false){
    ERROR_LOG << "You can only merge a HyperBinning with another HyperBinning" << std::endl;
    return;
  }
  
  const HyperBinning& otherHyperBinning = dynamic_cast<const HyperBinning&>(other);

  int nVolumes      =                   getNumHyperVolumes();
  int nVolumesOther = otherHyperBinning.getNumHyperVolumes();

  //this means every volume number in 'otherHyperBinning' needs to be increased by nVolumes.
  //This is important for linked bins and primary volume numbers!!
  
  reserveCapacity(nVolumes + nVolumesOther);

  for (int i = 0; i < nVolumesOther; i++){
    HyperVolume vol = otherHyperBinning.getHyperVolume(i);

    std::vector<int> linkedVolumes = otherHyperBinning.getLinkedHyperVolumes(i);

    for (unsigned int j = 0; j < linkedVolumes.size(); j++){
      linkedVolumes.at(j) += nVolumes;
    }

    addHyperVolume(vol, linkedVolumes);
  }
  
  int nPrimaryBinsOther = otherHyperBinning.getNumPrimaryVolumes();

  for (int i = 0; i < nPrimaryBinsOther; i++){
    int primaryVolumeNumber = otherHyperBinning.getPrimaryVolumeNumber(i);
    primaryVolumeNumber += nVolumes;
    addPrimaryVolumeNumber(primaryVolumeNumber);
  }

  _changed = true;

  //INFO_LOG << "Ending HyperBinning::mergeBinnings" << std::endl;


}





///Look at the tree that contains the HyperBinning and find the dimensionality
///
int HyperBinning::getHyperBinningDimFromTree(TTree* tree){

  if (tree == 0){
    ERROR_LOG << "Invalid tree in HyperBinning::getDimension(TTree* tree)" << std::endl;
    return 0;
  }  
  
  TString branchName = "lowCorner_0";
  int nDim = 0;

  while ( tree->GetListOfBranches()->FindObject(branchName) != 0 ){
    nDim++;
    branchName  = "lowCorner_";
    branchName += nDim;
  }
  
  if (nDim == 0){
    ERROR_LOG << "I cannot find any branches in the tree that indicate a HyperBinning is stored here" << std::endl;
    return 0;
  }

  return nDim;

}


///Create the branches in a TTree so that the HyperBinningMemRes
///can be saved.
void HyperBinning::createBranches(TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const{

  tree->Branch("binNumber", binNumber);
  tree->Branch("linkedBins", "vector<int>" ,linkedBins);
  for (int i = 0; i < getDimension(); i++) {
    TString lowCornerName  = "lowCorner_"; lowCornerName += i;
    TString highCornerName = "highCorner_"; highCornerName += i;
    tree->Branch(lowCornerName, lowCorner + i);
    tree->Branch(highCornerName, highCorner + i);
  }
  
}

///Save a single HyperVolume to a tree - this involves looping
///over every HyperCuboid in the HyperVolume.
void HyperBinning::saveHyperVolumeToTree(TTree* tree, double* lowCorner, double* highCorner, const HyperVolume& hyperVolume) const{

  for(int i = 0; i < hyperVolume.size(); i++){
    HyperCuboid hyperCuboid = hyperVolume.getHyperCuboid(i);
    HyperPoint lowCornerVect  = hyperCuboid.getLowCorner();
    HyperPoint highCornerVect = hyperCuboid.getHighCorner();
    for (int dim = 0; dim < getDimension(); dim++) lowCorner [dim] = lowCornerVect .at(dim);
    for (int dim = 0; dim < getDimension(); dim++) highCorner[dim] = highCornerVect.at(dim);
    tree->Fill();
  }

}

///Save the list of Primary Volume Numbers to the open (and in scope) TFile.
///
void HyperBinning::savePrimaryVolumeNumbers() const{

  TTree* tree = new TTree("PrimaryVolumeNumbers", "PrimaryVolumeNumbers");
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningMemRes::save()";
    return;
  }

  //Define branch addresses
  int volumeNumber = -1;

  tree->Branch("volumeNumber", &volumeNumber);
  
  int nPrimVols = getNumPrimaryVolumes();

  //Loop over each Primary Volume
  for(int i = 0; i < nPrimVols; i++ ){
    volumeNumber = getPrimaryVolumeNumber(i);
    tree->Fill();
  }
  
  //tree->Write();
  
}

///Save the HyperBinningMemRes to a TFile.
///
void HyperBinning::save(TString filename) const{

  TFile* file = new TFile(filename, "RECREATE");

  if (file == 0){
    ERROR_LOG << "Could not open TFile in HyperBinningMemRes::save(" << filename << ")";
    return;
  }

  save();

  //file->Write();
  file->Close();

}

///Save the HyperBinningMemRes to the open (and in scope) TFile.
///
void HyperBinning::save() const{
  
  savePrimaryVolumeNumbers();

  TTree* tree = new TTree("HyperBinning", "HyperBinning");
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningMemRes::save()";
    return;
  }

  //Define branch addresses
  int binNumber = -1;
  double* lowCorner = new double [getDimension()];
  double* highCorner = new double [getDimension()];
  std::vector<int>* linkedBins = new std::vector<int>();

  //Create branches and link them to branch addresses
  createBranches(tree, &binNumber, lowCorner, highCorner, &linkedBins);
  
  //Loop over each HyperVolume
  for(int bin = 0; bin < getNumHyperVolumes(); bin++ ){
    binNumber = bin;
    *linkedBins = getLinkedHyperVolumes(bin);
    //save all HyperCuboids in this HyperVolume to the TTree under the current bin number
    saveHyperVolumeToTree(tree, lowCorner, highCorner, getHyperVolume(bin));
  }
    
  delete lowCorner;
  delete highCorner;
  delete linkedBins;

}



std::vector<int> HyperBinning::getPrimaryVolumeNumbers() const{
  std::vector<int> primaryVolumeNumbers;
  int nPrimVols = getNumPrimaryVolumes();
  primaryVolumeNumbers.reserve(nPrimVols);

  for (int i = 0; i < nPrimVols; i++){
    primaryVolumeNumbers.push_back(getPrimaryVolumeNumber(i));
  }
  return primaryVolumeNumbers;
}




///Destructor
///
HyperBinning::~HyperBinning(){
  GOODBYE_LOG << "Goodbye from the HyperBinning() Constructor";
}




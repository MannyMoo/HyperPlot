#include "HyperBinningMemRes.h"


///The only constructor
HyperBinningMemRes::HyperBinningMemRes() 
{
  setBinningType("HyperBinningMemRes");
  WELCOME_LOG << "Hello from the HyperBinningMemRes() Constructor";
}

///Set the dimension of the HyperBinningMemRes. This can only be 
///called once, when it is known what dimesnion it is.
void HyperBinningMemRes::setDimension(int dim){
  
  if (getDimension() == 0){
    BinningBase::setDimension(dim);
    _averageBinWidth = HyperPoint ( getDimension(), 1.0);
    _minmax          = HyperCuboid( getDimension(), 0.0, 1.0 );
  }

}




/// This function will merge the two binnings. It assumes that the first
///
void HyperBinningMemRes::mergeBinnings( const BinningBase& other ){
  
  if (isSameBinningType(other) == false){
    ERROR_LOG << "You cannot merge a HyperBinningMemRes with a " << other.getBinningType() << std::endl;
    return;
  }
  
  const HyperBinningMemRes& otherHyperBinningMemRes = dynamic_cast<const HyperBinningMemRes&>(other);

  int nVolumes      = getNumHyperVolumes();
  int nVolumesOther = otherHyperBinningMemRes.getNumHyperVolumes();

  //this means every volume number in 'otherHyperBinningMemRes' needs to be increased by nVolumes.
  //This is important for linked bins and primary volume numbers!!

  for (int i = 0; i < nVolumesOther; i++){
    _hyperVolumes.push_back( otherHyperBinningMemRes._hyperVolumes.at(i) );
    std::vector<int> linkedVolumes = otherHyperBinningMemRes._linkedHyperVolumes.at(i);

    for (unsigned int j = 0; j < linkedVolumes.size(); j++){
      linkedVolumes.at(j) += nVolumes;
    }

    _linkedHyperVolumes.push_back(linkedVolumes);
  }
  
  int nPrimaryBinsOther = otherHyperBinningMemRes._primaryVolumeNumbers.size();

  for (int i = 0; i < nPrimaryBinsOther; i++){
    int primaryVolumeNumber = otherHyperBinningMemRes._primaryVolumeNumbers.at(i);
    primaryVolumeNumber += nVolumes;
    _primaryVolumeNumbers.push_back(primaryVolumeNumber);
  }

  _changed = true;

}

std::vector<int> HyperBinningMemRes::getLinkedHyperVolumes( int volumeNumber ) const{

  return _linkedHyperVolumes.at(volumeNumber);

}


BinningBase* HyperBinningMemRes::clone() const{

  return dynamic_cast<BinningBase*>(new HyperBinningMemRes(*this));

}



///Add a primary volume number
///
void HyperBinningMemRes::addPrimaryVolumeNumber(int volumeNumber){
  _primaryVolumeNumbers.push_back(volumeNumber);
}



///get the number of HyperVolumes
///
int HyperBinningMemRes::getNumHyperVolumes() const{
  return _hyperVolumes.size();
}  



int HyperBinningMemRes::getNumPrimaryVolumes  () const{
  return _primaryVolumeNumbers.size();
}

int HyperBinningMemRes::getPrimaryVolumeNumber(int i) const{
  return _primaryVolumeNumbers.at(i);
}

///Add a HyperVolume to the HyperBinningMemRes and add a set of empty
///HyperVolume links.
bool HyperBinningMemRes::addHyperVolume(const HyperVolume& hyperVolume, std::vector<int> linkedVolumes){
  
  //If this is the first volume that has been added, use it to set the dimension
  if (_hyperVolumes.size() == 0){
    setDimension( hyperVolume.getDimension() );
  }

  if (hyperVolume.getDimension() == getDimension()) {
    _hyperVolumes.push_back(hyperVolume); 
    _linkedHyperVolumes.push_back(linkedVolumes);
    _changed = true;
    return true;
  }

  ERROR_LOG << "This HyperVolume has the wrong dimensionality for this HyperBinningMemRes";
  return false;

}


///Create the branches in a TTree so that the HyperBinningMemRes
///can be saved.
void HyperBinningMemRes::createBranches(TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const{

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
void HyperBinningMemRes::saveHyperVolumeToTree(TTree* tree, double* lowCorner, double* highCorner, const HyperVolume& hyperVolume) const{

  for(int i = 0; i < hyperVolume.size(); i++){
    HyperCuboid hyperCuboid = hyperVolume.getHyperCuboid(i);
    HyperPoint lowCornerVect  = hyperCuboid.getLowCorner();
    HyperPoint highCornerVect = hyperCuboid.getHighCorner();
    for (int dim = 0; dim < getDimension(); dim++) lowCorner [dim] = lowCornerVect .at(dim);
    for (int dim = 0; dim < getDimension(); dim++) highCorner[dim] = highCornerVect.at(dim);
    tree->Fill();
  }

}

///Save the HyperBinningMemRes to a TFile.
///
void HyperBinningMemRes::save(TString filename) const{

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
void HyperBinningMemRes::save() const{
  
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
  for(unsigned int bin = 0; bin < _hyperVolumes.size(); bin++ ){
    binNumber = bin;
    *linkedBins = _linkedHyperVolumes.at(bin);
    //save all HyperCuboids in this HyperVolume to the TTree under the current bin number
    saveHyperVolumeToTree(tree, lowCorner, highCorner, _hyperVolumes.at(bin));
  }
  
  //tree->Write();
  
  delete lowCorner;
  delete highCorner;
  delete linkedBins;

}

///Save the list of Primary Volume Numbers to the open (and in scope) TFile.
///
void HyperBinningMemRes::savePrimaryVolumeNumbers() const{

  TTree* tree = new TTree("PrimaryVolumeNumbers", "PrimaryVolumeNumbers");
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningMemRes::save()";
    return;
  }

  //Define branch addresses
  int volumeNumber = -1;

  tree->Branch("volumeNumber", &volumeNumber);

  //Loop over each Primary Volume
  for(unsigned int i = 0; i < _primaryVolumeNumbers.size(); i++ ){
    volumeNumber = _primaryVolumeNumbers.at(i);
    tree->Fill();
  }
  
  //tree->Write();
  
}

///Save the list of Primary Volume Numbers to the open (and in scope) TFile.
///
void HyperBinningMemRes::loadPrimaryVolumeNumbers(TFile* file) {

  TTree* tree = dynamic_cast<TTree*>( file->Get("PrimaryVolumeNumbers") );
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningMemRes::loadPrimaryVolumeNumbers()";
    return;
  }

  //Define branch addresses
  int volumeNumber = -1;

  tree->SetBranchAddress("volumeNumber", &volumeNumber);

  //Loop over each Primary Volume
  for(int i = 0; i < tree->GetEntries(); i++ ){
    tree->GetEntry(i);
    _primaryVolumeNumbers.push_back(volumeNumber);
  }
  
}

///Set branch addresses for loading HyperBinningMemRes
///from a file.
void HyperBinningMemRes::setBranchAddresses(TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const{

  tree->SetBranchAddress("binNumber", binNumber);
  tree->SetBranchAddress("linkedBins", linkedBins);
  for (int i = 0; i < getDimension(); i++) {
    TString lowCornerName  = "lowCorner_"; lowCornerName += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;
    tree->SetBranchAddress(lowCornerName, lowCorner + i);
    tree->SetBranchAddress(highCornerName, highCorner + i);
  }

}


HyperVolume HyperBinningMemRes::getHyperVolume(int volumeNumber) const{
  return _hyperVolumes.at(volumeNumber);
}


///Look at the tree that contains the HyperBinningMemRes and find the dimensionality
///
std::vector<int> HyperBinningMemRes::getPrimaryVolumeNumbers() const{
  return _primaryVolumeNumbers;
}


///Load HyperBinningMemRes from a file
///
void HyperBinningMemRes::load(TString filename, TString option){
  
  if (option != "READ"){
    INFO_LOG << "For a memory resident HyperBinning you should always use the READ option. Setting to READ" << std::endl;
    option = "READ";
  }

  TFile* file = new TFile(filename, "READ");

  if (file == 0){
    ERROR_LOG << "Could not open TFile in HyperBinningMemRes::load(" << filename << ")";
    return;
  }

  loadPrimaryVolumeNumbers(file);

  TTree* tree = (TTree*)file->Get("HyperBinning");

  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningMemRes::load()";
    return;
  }
  
  //Figure out how many dimensions there are from the tree
  setDimension( getHyperBinningDimFromTree(tree) );

  //Create branch addresses and link them to TTree
  int binNumber = -1;
  double* lowCorner  = new double [getDimension()];
  double* highCorner = new double [getDimension()];  
  std::vector<int>* linkedBins = new std::vector<int>();

  setBranchAddresses(tree, &binNumber, lowCorner, highCorner, &linkedBins);
  

  //Loop over the TTree and fill the HyperBinningMemRes
  int nEntries = tree->GetEntries();

  int currentBinNumber = -1;
  HyperVolume* currentHyperVolume = new HyperVolume(getDimension());
  std::vector<int> currentLinkedVolumes;

  //_verbose = true;

  for(int ent = 0; ent < nEntries; ent++){
    tree->GetEntry(ent);
    
    //If the bin number hasn't changed, need to
    //add HyperCube to previous HyperVolume
    if(ent == 0 || currentBinNumber == binNumber){
      currentBinNumber = binNumber;
      currentLinkedVolumes = *linkedBins;
      HyperPoint lowCornerVect (getDimension());
      HyperPoint highCornerVect(getDimension());
      for (int i = 0; i < getDimension(); i++){
        lowCornerVect .at(i) = lowCorner [i];
        highCornerVect.at(i) = highCorner[i];
      }
      VERBOSE_LOG << "Adding cuboid to volume";
      currentHyperVolume->addHyperCuboid(lowCornerVect, highCornerVect);
    }

    //if the bin number has changed, need to add
    //previous HyperVolume to HyperBinningMemRes
    //and start a new HyperVolume
    else {
      VERBOSE_LOG << "Adding volume to binning";
      this->addHyperVolume(*currentHyperVolume, currentLinkedVolumes);
      delete currentHyperVolume;
      currentHyperVolume = new HyperVolume(getDimension());
      currentLinkedVolumes = *linkedBins;
      currentBinNumber = binNumber;
      HyperPoint lowCornerVect (getDimension());
      HyperPoint highCornerVect(getDimension());
      for (int i = 0; i < getDimension(); i++){
        lowCornerVect .at(i) = lowCorner [i];
        highCornerVect.at(i) = highCorner[i];
      }
      VERBOSE_LOG << "Adding cuboid to volume";
      currentHyperVolume->addHyperCuboid(lowCornerVect, highCornerVect);
      VERBOSE_LOG << "Adding linked volumes";
      
    }

    //if it's the final iteration, need to add the final HyperVolume
    if (ent == nEntries - 1) {
      this->addHyperVolume(*currentHyperVolume);
    }
  }
  
  VERBOSE_LOG << "Binning loaded";

  //_verbose = false;

  delete currentHyperVolume;
  delete lowCorner;
  delete highCorner;
  delete linkedBins;

  _changed = true;
  
  file->Close();

}

///Destructor
///
HyperBinningMemRes::~HyperBinningMemRes(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMemRes() Constructor";
}




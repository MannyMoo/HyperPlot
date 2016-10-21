#include "HyperBinningDiskRes.h"


///The only constructor
HyperBinningDiskRes::HyperBinningDiskRes() :
  _file(0),
  _writeable(false),
  _tree(0),
  _cuboid(0),
  _linkedBins(new std::vector<int>()),
  _volumeNumber(-1),
  _currentEntry(-1)
{
  setBinningType("HyperBinningDiskRes");
  WELCOME_LOG << "Hello from the HyperBinningDiskRes() Constructor";
}

///Set the dimension of the HyperBinningDiskRes. This can only be 
///called once, when it is known what dimesnion it is.
void HyperBinningDiskRes::setDimension(int dim){
  
  if (getDimension() == 0){
    HyperBinning::setDimension(dim);
    _cuboid = HyperCuboid(dim);
  }

}


void HyperBinningDiskRes::getEntry(int volumeNumber) const{
  if (_tree == 0){
    ERROR_LOG << "HyperBinningDiskRes::getEntry(int volumeNumber), tree doesn't exist" << std::endl; 
    return;
  }

  if ( _currentEntry != volumeNumber ){
    _currentEntry = volumeNumber;
    _tree->GetEntry(volumeNumber);
  }
}

HyperVolume HyperBinningDiskRes::getHyperVolume(int volumeNumber) const{
  getEntry(volumeNumber);
  return HyperVolume(_cuboid);
} 

std::vector<int> HyperBinningDiskRes::getLinkedHyperVolumes( int volumeNumber ) const{
  getEntry(volumeNumber);
  return *_linkedBins;
}




/// This function will merge the two binnings. It assumes that the first
///
void HyperBinningDiskRes::mergeBinnings( const BinningBase& other ){
  
  if ( other.getBinningType().Contains("HyperBinning") == false){
    ERROR_LOG << "You can only merge a HyperBinningDiskRes with another HyperBinning" << std::endl;
    return;
  }
  
  const HyperBinning& otherHyperBinning = dynamic_cast<const HyperBinning&>(other);

  int nVolumes      = getNumHyperVolumes();
  int nVolumesOther = otherHyperBinning.getNumHyperVolumes();

  //this means every volume number in 'otherHyperBinning' needs to be increased by nVolumes.
  //This is important for linked bins and primary volume numbers!!

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
    int primaryVolumeNumber = getPrimaryVolumeNumber(i);
    primaryVolumeNumber += nVolumes;
    addPrimaryVolumeNumber(primaryVolumeNumber);
  }

  _changed = true;

}




BinningBase* HyperBinningDiskRes::clone() const{

  return dynamic_cast<BinningBase*>(new HyperBinningDiskRes(*this));

}


///Add a primary volume number
///
void HyperBinningDiskRes::addPrimaryVolumeNumber(int volumeNumber){

  if (_writeable == false){
    ERROR_LOG << "You cannot wirte to this Disk Resident HyperBinning - you must have opened in READ mode" << std::endl;
    return;
  }

  _primVolNum = volumeNumber;
  _treePrimVol->Fill();
}



///get the number of HyperVolumes
///
int HyperBinningDiskRes::getNumHyperVolumes() const{
  if (_tree == 0){
    return 0;
  }
  return _tree->GetEntries();
}  





///Save the HyperBinningDiskRes to a TFile.
///
void HyperBinningDiskRes::save(TString filename) const{
  ERROR_LOG << "Cannot save a Disk Resident HyperBinning to " << filename;
  ERROR_LOG << " becasue it already exists on disk!" << std::endl;
}

///Save the HyperBinningDiskRes to the open (and in scope) TFile.
///
void HyperBinningDiskRes::save() const{
  

}

bool HyperBinningDiskRes::addHyperVolume(const HyperVolume& hyperVolume, std::vector<int> linkedVolumes){
  
  if (_writeable == false){
    ERROR_LOG << "You cannot wirte to this Disk Resident HyperBinning - you must have opened in READ mode" << std::endl;
    return false;
  }
  if (getDimension() == 0){
    setDimension(hyperVolume.getDimension());
    createHyperBinningTree();
  }
  
  if (_tree == 0){
    ERROR_LOG << "HyperBinningDiskRes::addHyperVolume - there isn't a tree to write to" << std::endl;
    return false;
  }

  int nVolumes = getNumHyperVolumes();
  for (int i = 0; i < hyperVolume.size(); i++){
    *_linkedBins  = linkedVolumes;
    _cuboid       = hyperVolume.at(i);
    _volumeNumber = nVolumes;
    _tree->Fill();
  }
  
  return true;
}




int HyperBinningDiskRes::getNumPrimaryVolumes  () const{
  if (_treePrimVol == 0) return 0;
  return _treePrimVol->GetEntries();
}

int HyperBinningDiskRes::getPrimaryVolumeNumber(int i) const{
  if (_treePrimVol == 0) {
    ERROR_LOG << "HyperBinningDiskRes::getPrimaryVolumeNumber(i). No primary volume tree" << std::endl;
  }
  _treePrimVol->GetEntry(i);
  return _primVolNum;
}


void HyperBinningDiskRes::loadHyperBinningTree(){

  _tree = dynamic_cast<TTree*>(_file->Get("HyperBinning"));

  if (_tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningDiskRes::load()";
    return;
  }
  
  //Figure out how many dimensions there are from the tree
  setDimension( getHyperBinningDimFromTree(_tree) );
  
  _tree->SetBranchAddress("binNumber" , &_volumeNumber);
  _tree->SetBranchAddress("linkedBins", &_linkedBins   );

  for (int i = 0; i < getDimension(); i++){
    TString lowCornerName  = "lowCorner_" ; lowCornerName  += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;   
    _tree->SetBranchAddress(lowCornerName , &_cuboid.getLowCorner ().at(i) );
    _tree->SetBranchAddress(highCornerName, &_cuboid.getHighCorner().at(i) );
  }
  
  _currentEntry = -1;
  _changed = true;

}



void HyperBinningDiskRes::loadPrimaryVolumeTree(){

  _treePrimVol = dynamic_cast<TTree*>( _file->Get("PrimaryVolumeNumbers") );
  
  if (_treePrimVol == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningDiskRes::loadPrimaryVolumeNumbers()";
    return;
  }

  _treePrimVol->SetBranchAddress("volumeNumber", &_primVolNum);

}

void HyperBinningDiskRes::createHyperBinningTree(){
  
  INFO_LOG << "Creating a tree to store DiskRes HyperBinning" << std::endl;

  _file->cd();
  _tree = new TTree("HyperBinning", "HyperBinning");

  //Figure out how many dimensions there are from the tree
  int dim = getDimension();
  if (dim == 0){
    ERROR_LOG << "The dimesion has not yet been set, so I cannot createHyperBinningTree!!" << std::endl;
  }

  _tree->Branch("binNumber" , &_volumeNumber);
  _tree->Branch("linkedBins", &_linkedBins   );

  for (int i = 0; i < dim; i++){
    TString lowCornerName  = "lowCorner_" ; lowCornerName  += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;   
    _tree->Branch(lowCornerName , &_cuboid.getLowCorner ().at(i) );
    _tree->Branch(highCornerName, &_cuboid.getHighCorner().at(i) );
  }
  
  _currentEntry = -1;
  _changed = true;

}

void HyperBinningDiskRes::createPrimaryVolumeTree(){

  _file->cd();
  _treePrimVol = new TTree("PrimaryVolumeNumbers", "PrimaryVolumeNumbers");
  _treePrimVol->Branch("volumeNumber", &_primVolNum);

}

///Load HyperBinningDiskRes from a file
///
void HyperBinningDiskRes::load(TString filename, TString option){
  
  if (option == "UPDATE" || option == "RECREATE") {_writeable = true;}

  _file = new TFile(filename, option);
  
  if (_file == 0){
    ERROR_LOG << "Could not open TFile in HyperBinningDiskRes::load(" << filename << ")";
    return;
  }
  
  if (option == "UPDATE" || option == "READ"){
    loadHyperBinningTree   ();
    loadPrimaryVolumeTree  ();
  }
  else if (option == "RECREATE"){
    //createHyperBinningTree   ();
    createPrimaryVolumeTree  ();    
  }
  else{
    ERROR_LOG << "There are only three load options (READ, RECREATE, UPDATE) - You have selected " << option << std::endl;
    return;
  }

  INFO_LOG << "Sucessfully attached Disk Resident HyperBinning to file " << filename << std::endl;

}

bool HyperBinningDiskRes::isDiskResident() const{
  return true;
}
TString HyperBinningDiskRes::filename() const{
  return _file->GetName();
}


///Destructor
///
HyperBinningDiskRes::~HyperBinningDiskRes(){
  GOODBYE_LOG << "Goodbye from the HyperBinningDiskRes() Constructor";
  if (_file != 0) {
    INFO_LOG << "Closing HyperBinningDiskRes " << _file->GetName() << std::endl;
    _file->cd();
    if (_writeable == true){
      _tree->Write();
      _treePrimVol->Write();
    }
    _file->Close();
  }
}




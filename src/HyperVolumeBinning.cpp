#include "HyperVolumeBinning.h"

///The only constructor, just takes the dimensionality of
///the binning
HyperVolumeBinning::HyperVolumeBinning(int dimension) :
  _dimension(dimension),
  _changed(true),
  _averageBinWidth(dimension, 1.0),
  _minmax( HyperPoint(dimension, 0.0), HyperPoint(dimension, 1.0) ),
  _names(dimension, "")
{
  WELCOME_LOG << "Hello from the HyperVolumeBinning() Constructor";
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
int HyperVolumeBinning::getBinNum(const HyperPoint& coords) const{
  
  if (_changed == true) updateCash();

  //First check if the HyperPoint is in the HyperCuboid _minmax that
  //surrounds all the bins.

  if (_minmax.inVolume(coords) == 0) return -1;
  
  
  if ( _primaryVolumeNumbers.size() == 0){

    //loop over all the bins until one is found that contains the event. If a bin
    // is found that contains the event, check if it has any linked bins. Linked bins
    // aren't real bins... just used to speed up sorting later.

    int volumeNumber = -1;
  
    for (unsigned int i = 0; i < _hyperVolumes.size(); i++){
      bool inVol = _hyperVolumes.at(i).inVolume(coords);
      if (inVol == 1) { volumeNumber = i; break; }
    }
     
    if (volumeNumber == -1) return -1;
  
    if ( _linkedHyperVolumes.at(volumeNumber).size() > 0 ) volumeNumber = followBinLinks(coords, volumeNumber);
  
    return _binNum.at(volumeNumber);
  }


  int primaryVolumeNumber = -1;

  for (unsigned int i = 0; i < _primaryVolumeNumbers.size(); i++){
    int thisVolNum = _primaryVolumeNumbers.at(i);
    bool inVol = _hyperVolumes.at(thisVolNum).inVolume(coords);
    if (inVol == 1) { primaryVolumeNumber = thisVolNum; break; }
  }
  
  int volumeNumber = -1;

  if ( _linkedHyperVolumes.at(primaryVolumeNumber).size() > 0 ) {
    volumeNumber = followBinLinks(coords, primaryVolumeNumber);
  }
  else{
    ERROR_LOG << "This primary volume has NO links. Not what I expect!!" << std::endl;
    return _binNum.at(primaryVolumeNumber);
  }
  
  return _binNum.at(volumeNumber);

}


/// This function will merge the two binnings. It assumes that the first
///
void HyperVolumeBinning::mergeBinnings( const HyperVolumeBinning& other ){

  int nVolumes      = getNumHyperVolumes();
  int nVolumesOther = other.getNumHyperVolumes();

  //this means every volume number in 'other' needs to be increased by nVolumes.
  //This is important for linked bins and primary volume numbers!!

  for (int i = 0; i < nVolumesOther; i++){
    _hyperVolumes.push_back( other._hyperVolumes.at(i) );
    std::vector<int> linkedVolumes = other._linkedHyperVolumes.at(i);

    for (unsigned int j = 0; j < linkedVolumes.size(); j++){
      linkedVolumes.at(j) += nVolumes;
    }

    _linkedHyperVolumes.push_back(linkedVolumes);
  }
  
  int nPrimaryBinsOther = other._primaryVolumeNumbers.size();

  for (int i = 0; i < nPrimaryBinsOther; i++){
    int primaryVolumeNumber = other._primaryVolumeNumbers.at(i);
    primaryVolumeNumber += nVolumes;
    _primaryVolumeNumbers.push_back(primaryVolumeNumber);
  }

  _changed = true;

}



/// \todo This doesn't really work yet
///
void HyperVolumeBinning::drawBin(int dim, int bin, RootPlotter2D& plotter) const{

  double min = getBinHyperVolume(bin).getHyperCuboid(0).getLowCorner() .at(dim); 
  double max = getBinHyperVolume(bin).getHyperCuboid(0).getHighCorner().at(dim); 

  TLine* line = new TLine( bin , min, bin, max );

  double nBins  = this->getNumBins();
  double lineWidth = (50.0/nBins)*8.0;
  
  if ( lineWidth > 10.0 ) lineWidth = 10.0;

  line->SetLineWidth(lineWidth);

  plotter.addObject(line);

}

/// \todo This doesn't really work yet
///
void HyperVolumeBinning::drawBinning(int dim, TString name, TPad* pad, double upperMargin, double lowerMargin) const{
  
  int nBins  = this->getNumBins();
  double min = this->getMin(dim);
  double max = this->getMax(dim);
  
  int extra = ceil(nBins*0.05);
  double extraD = (max - min)*0.05;
  
  TString uniqueName = "histBinningDim"; uniqueName += dim;
  TString axisTitle = "dim "; axisTitle += dim;

  TH2D* hist = new TH2D(uniqueName, uniqueName, nBins + extra*2.0, -extra, (double)nBins - 1 + extra, 10, min - extraD, max + extraD);
  
  hist->GetYaxis()->CenterTitle();
  //hist->GetXaxis()->CenterTitle();

  hist->GetXaxis()->SetTitle("Bin Number");
  hist->GetYaxis()->SetTitle("dim");

  RootPlotter2D plotter(hist);
  
  plotter.scaleTextSize(2.0);

  plotter.setBMargin(lowerMargin);
  plotter.setTMargin(upperMargin);
  plotter.setYAxisTitleOffset(0.4);

  for (int i = 0; i < nBins; i++) drawBin(dim, i, plotter);

  double usualMargins = 0.78;
  double margins = 1.0 - upperMargin - lowerMargin;

  double scale = margins/usualMargins;
  scale *= 0.9;

  plotter.plot(name, "", pad, scale );  

  delete hist;  

}

/// \todo This doesn't really work yet
///
void HyperVolumeBinning::drawBinning(TString name, TPad* pad) const{

  int nDim = this->getDimension();

  double indivPadHeight = 80.0; 
  double extraPadHegiht = 30.0; //for bottom and top pad with axis
  
  double marginFrac = extraPadHegiht/(extraPadHegiht+indivPadHeight);

  double canvasheight = indivPadHeight*nDim + 2.0*extraPadHegiht;
  
  
  TPad* canvas = 0;

  if (pad == 0) canvas = new TCanvas("canvas", "canvas",400, canvasheight);
  else canvas = pad;



  double extraPadFrac = extraPadHegiht / canvasheight;
  double indivPadFrac = indivPadHeight / canvasheight;

 
  
  double lower = 0.01;
  double upper = extraPadFrac + indivPadFrac;

  for(int i = nDim - 1; i >= 0; i--){

    TString uniquePadName = name;  uniquePadName += i;
    canvas->cd();
    TPad* pad = new TPad(uniquePadName, uniquePadName, 0.01, lower, 0.99, upper);
    
    lower = upper;
    upper += indivPadFrac;
    if (i == 0) upper = 0.99;

    double upMargin = 0.01; 
    double loMargin = 0.01; 
    if (i == 0       ) upMargin = marginFrac;
    if (i == nDim - 1) loMargin = marginFrac;

    drawBinning(i, "", pad ,upMargin, loMargin);
    canvas->cd();
    pad->Draw();
  }
  
  //if (name == "") canvas->Draw();
  if (name != "") canvas->Print(name + Plotter::s_imageformat);

}

///Add a primary volume number
///
void HyperVolumeBinning::addPrimaryVolumeNumber(int volumeNumber){
  _primaryVolumeNumbers.push_back(volumeNumber);
}


///Used to follow the bin hierarchy. Give it a HyperPoint, and the number of a 
///HyperVolume (that has links) that the HyperPoint falls into. 
///
///
int HyperVolumeBinning::followBinLinks(const HyperPoint& coords, int motherVolumeNumber) const{
  
  //find the linked volumes
  const std::vector<int>& linkedVolumes = _linkedHyperVolumes.at(motherVolumeNumber);
  
  int volumeNumber = -1;
  
  //see if the coords falls into any of the linked volumes (it should if there are no bugs)
  for (unsigned i = 0; i < linkedVolumes.size(); i++){
    int daughBinNum = linkedVolumes.at(i);
    bool inVol = _hyperVolumes.at(daughBinNum).inVolume(coords);
    if (inVol == 1) { volumeNumber = daughBinNum; break; }
  }
  
  if (volumeNumber == -1) {
    ERROR_LOG << "The trail of linked bins has gone cold!";
    return -1;
  }
  
  //now have volumeNumber which contains the next bin in the hierarchy.
  // if this is linked to more bins, keep following the trail!
  if ( _linkedHyperVolumes.at(volumeNumber).size() > 0 ) volumeNumber = followBinLinks(coords, volumeNumber);
  
  //if not, we have made it to the end. Return the volume number!
  return volumeNumber;

}


///Look at the distance between the given HyperPoint and the 'AverageCenter' (or CenterOfMass)
///of all the HyperVolume bins. Return the n that are closest.
///
std::vector<int> HyperVolumeBinning::findNNearestBins( const HyperPoint& point, int n ) const{

  int        * arrHyperVol = new int    [getNumBins()];
  double     * arrDist     = new double [getNumBins()];  
  int        * arrIndex    = new int    [getNumBins()]; 

  int count = 0;

  for (unsigned int i = 0; i < _hyperVolumes.size(); i++){
    if ( getBinNum(i) == -1 ) continue;
    HyperPoint center = _hyperVolumes.at(i).getAverageCenter();
    double dist = center.distanceTo(point);
    arrHyperVol[count] = i;
    arrDist    [count] = dist;
    count++;
  }

  TMath::Sort( getNumBins(), arrDist, arrIndex, false);

  std::vector<int> nearest;
  for (int i = 0; i < n; i++){
    nearest.push_back( arrHyperVol[ arrIndex[i] ] );
    //std::cout << "   " << arrHyperVol[ arrIndex[i] ] << "     " << getBinNum( arrHyperVol[ arrIndex[i] ] ) << std::endl;
  }

  delete arrIndex;
  delete arrDist;

  return nearest;

}


///Imagine a line that starts from the given HyperPoint, and travels in the dimension given.
///This function finds all HyperVolume bins that intersect this line, and then returns the 
///bin numbers. The bin numbers are ordered along the direction of the line.
std::vector<int> HyperVolumeBinning::findOrderedBinsOnLine( const HyperPoint& point, int dim ) const{

  std::vector<int> volumesOnLine;
  std::vector<double> binCenters;

  for (unsigned int i = 0; i < _hyperVolumes.size(); i++){
    if ( _binNum.at(i) == -1 ) continue;
    if ( isLineInVolume( _hyperVolumes.at(i), point, dim ) == true){
      volumesOnLine.push_back( i );
      binCenters   .push_back( _hyperVolumes.at(i).getAverageCenter().at(dim) );
    }
  }
  
  int        * arrIndex = new int    [binCenters.size()];
  double     * arrCen   = new double [binCenters.size()];

  for (unsigned i = 0; i < binCenters.size(); i++){
    arrIndex[i] = i;
    arrCen  [i] = binCenters.at(i);
  }
  
  TMath::Sort( (int)binCenters.size(), arrCen, arrIndex, false);
  
  std::vector<int> sortedVolumesOnLine;

  for (unsigned i = 0; i < binCenters.size(); i++){
    sortedVolumesOnLine.push_back( volumesOnLine.at( arrIndex[i] ) );
  }
  
  delete [] arrIndex;
  delete [] arrCen;

  return sortedVolumesOnLine;

}

///Check if a line that starts from the given HyperPoint, and travels in the dimension given,
///goes through the HyperVolume given.
///
bool HyperVolumeBinning::isLineInVolume( const HyperVolume& volume, const HyperPoint& point, int dim ) const{
  
  int nCubes = volume.getHyperCuboids().size();

  for (int i = 0; i < nCubes; i++){
    bool val = isLineInCuboid( volume.getHyperCuboid(i), point, dim );
    if (val == true) return true;
  }
  
  return false;

}

///Check if a line that starts from the given HyperPoint, and travels in the dimension given,
///goes through the HyperCuboid given.
///
bool HyperVolumeBinning::isLineInCuboid( const HyperCuboid& cuboid, const HyperPoint& point, int dim ) const{

  for (int i = 0; i < getDimension(); i++){
    if ( i == dim ) continue;
    HyperPoint lowCorner  = cuboid.getLowCorner ();
    HyperPoint highCorner = cuboid.getHighCorner();
    if ( lowCorner.at(i) > point.at(i) || highCorner.at(i) <= point.at(i) ) return false;
  }
  
  return true;

}


///Get number of bins (this is NOT the number of
///HyperVolumes!!! - see the class description for more details)
int HyperVolumeBinning::getNumBins() const{
  
  //std::cout << "getNumBins" << std::endl;
  if (_changed == true) updateCash();
  
  return _hyperVolumeNumFromBinNum.size();

}

///get the number of HyperVolumes
///
int HyperVolumeBinning::getNumHyperVolumes() const{
  return _hyperVolumes.size();
}  

///Get the HyperVolume assosiated with bin number i. This just uses 
///the _hyperVolumeNumFromBinNum variable to find the HyperVolume number
///from the bin number, then returns that HyperVolume.
const HyperVolume& HyperVolumeBinning::getBinHyperVolume(int binNumber) const{
  
  //std::cout << "getBinHyperVolume" << std::endl;
  if (_changed == true) updateCash();
  return _hyperVolumes.at( _hyperVolumeNumFromBinNum.at(binNumber) );

}

///Get the bin number assosiated with a given HyperVolume number. 
///If this returns -1, it means that the HyperVolume in question
///is not a bin, but part of the binning hierarchy.
int HyperVolumeBinning::getBinNum(int volumeNumber) const{
  if (_changed == true) updateCash();
  return _binNum.at(volumeNumber);
}

/// get the HyperVolume Number from the bin number
///
int HyperVolumeBinning::getHyperVolumeNumber(int binNumber) const{
  if (_changed == true) updateCash(); 
  return _hyperVolumeNumFromBinNum.at(binNumber);
}

///Update the cash which includes the  mutable member variables
///_binNum, _hyperVolumeNumFromBinNum, _averageBinWidth,
/// and _minmax.
void HyperVolumeBinning::updateCash() const{

  //note the ordering of theses is important...
  //possilbe infinite loops

  _changed = false;

  updateBinNumbering();
  updateAverageBinWidth();
  updateMinMax();

}

///Update the member variables _binNum and _hyperVolumeNumFromBinNum.
///Will usually be called from updateCash()
void HyperVolumeBinning::updateBinNumbering() const{
  
  //first fill all the bin numbers with -1
  int nVolumes = getNumHyperVolumes();
  _binNum = std::vector<int>(nVolumes, -1);

  //if a HyperVolume has any linked HyperVolumes,
  //then set its bin number to count.
  int count = 0;
  for (int i = 0; i < getNumHyperVolumes(); i++){
    if ( _linkedHyperVolumes.at(i).size() == 0 ) {
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

///update the _averageBinWidth HyperPoint. 
///Will usually be called from updateCash()
void HyperVolumeBinning::updateAverageBinWidth() const{
  
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
void HyperVolumeBinning::updateMinMax() const{
  
  int dim = getDimension();

  HyperPoint min(dim);
  HyperPoint max(dim);
  
  for (int d = 0; d < dim; d++){
    min.at(d) = _hyperVolumes.at(0).getMin(d);
    max.at(d) = _hyperVolumes.at(0).getMax(d);
  }

  for(unsigned int i = 1; i < _hyperVolumes.size(); i++){
    const HyperVolume& thisVol = _hyperVolumes.at(i);
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
HyperPoint HyperVolumeBinning::getAverageBinWidth() const{
  //std::cout << "getAverageBinWidth" << std::endl;
  if (_changed == true) updateCash();
  return _averageBinWidth;  

}

///return the minimum value for a selected dimension.
///This value is cashed for speed - when the binning changes the cashe will
///automatically be updated.
double HyperVolumeBinning::getMin(int dimension) const{

  if (_changed == true) updateCash();
  return _minmax.getLowCorner().at(dimension);

}
///return the maximum value for a selected dimension.
///This value is cashed for speed - when the binning changes the cashe will
///automatically be updated.
double HyperVolumeBinning::getMax(int dimension) const{

  if (_changed == true) updateCash();
  return _minmax.getHighCorner().at(dimension);

}


///Add a HyperVolume to the HyperVolumeBinning and add a set of empty
///HyperVolume links.
bool HyperVolumeBinning::addHyperVolume(const HyperVolume& hyperVolume){

  if (hyperVolume.getDimension() == _dimension) {
    _hyperVolumes.push_back(hyperVolume); 
    _linkedHyperVolumes.push_back(std::vector<int>(0, 0.0));
    _changed = true;
    return true;
  }

  ERROR_LOG << "This HyperVolume has the wrong dimensionality for this HyperVolumeBinning";
  return false;

}

///Add a link from one HyperVolume to another HyperVolume. 
///Used for the hierarchy of bins discussed in the class description.
bool HyperVolumeBinning::addHyperVolumeLink(int volumeNum, int linkedVolumeNum){

  if ( volumeNum < getNumHyperVolumes() && volumeNum >= 0) {
    _linkedHyperVolumes.at(volumeNum).push_back(linkedVolumeNum);
    _changed = true;
    return true;
  }

  ERROR_LOG << "The HyperVolume you are trying to link to does not exist";
  return false;
}

///Create the branches in a TTree so that the HyperVolumeBinning
///can be saved.
void HyperVolumeBinning::createBranches(TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const{

  tree->Branch("binNumber", binNumber);
  tree->Branch("linkedBins", "vector<int>" ,linkedBins);
  for (int i = 0; i < _dimension; i++) {
    TString lowCornerName  = "lowCorner_"; lowCornerName += i;
    TString highCornerName = "highCorner_"; highCornerName += i;
    tree->Branch(lowCornerName, lowCorner + i);
    tree->Branch(highCornerName, highCorner + i);
  }
  
}

///Save a single HyperVolume to a tree - this involves looping
///over every HyperCuboid in the HyperVolume.
void HyperVolumeBinning::saveHyperVolumeToTree(TTree* tree, double* lowCorner, double* highCorner, const HyperVolume& hyperVolume) const{

  for(int i = 0; i < hyperVolume.size(); i++){
    HyperCuboid hyperCuboid = hyperVolume.getHyperCuboid(i);
    HyperPoint lowCornerVect  = hyperCuboid.getLowCorner();
    HyperPoint highCornerVect = hyperCuboid.getHighCorner();
    for (int dim = 0; dim < _dimension; dim++) lowCorner [dim] = lowCornerVect .at(dim);
    for (int dim = 0; dim < _dimension; dim++) highCorner[dim] = highCornerVect.at(dim);
    tree->Fill();
  }

}

///Save the HyperVolumeBinning to a TFile.
///
void HyperVolumeBinning::save(TString filename) const{

  TFile* file = new TFile(filename, "RECREATE");

  if (file == 0){
    ERROR_LOG << "Could not open TFile in HyperVolumeBinning::save(" << filename << ")";
    return;
  }

  save();

  file->Write();
  file->Close();

}

///Save the HyperVolumeBinning to the open (and in scope) TFile.
///
void HyperVolumeBinning::save() const{
  
  savePrimaryVolumeNumbers();

  TTree* tree = new TTree("HyperVolumeBinning", "HyperVolumeBinning");
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperVolumeBinning::save()";
    return;
  }

  //Define branch addresses
  int binNumber = -1;
  double* lowCorner = new double [_dimension];
  double* highCorner = new double [_dimension];
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
  
  tree->Write();
  
  delete lowCorner;
  delete highCorner;
  delete linkedBins;

}

///Save the list of Primary Volume Numbers to the open (and in scope) TFile.
///
void HyperVolumeBinning::savePrimaryVolumeNumbers() const{

  TTree* tree = new TTree("PrimaryVolumeNumbers", "PrimaryVolumeNumbers");
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperVolumeBinning::save()";
    return;
  }

  //Define branch addresses
  int volumeNumber = -1;

  tree->Branch("volumeNumber", &volumeNumber);

  //Loop over each Primary Volume
  for(unsigned int i = 0; i < _primaryVolumeNumbers.size(); i++ ){
    volumeNumber = _primaryVolumeNumbers.at(i);
    tree->Write();
  }
  
  tree->Write();
  
}

///Save the list of Primary Volume Numbers to the open (and in scope) TFile.
///
void HyperVolumeBinning::loadPrimaryVolumeNumbers(TFile* file) {

  TTree* tree = dynamic_cast<TTree*>( file->Get("PrimaryVolumeNumbers") );
  
  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperVolumeBinning::loadPrimaryVolumeNumbers()";
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

///Set branch addresses for loading HyperVolumeBinning
///from a file.
void HyperVolumeBinning::setBranchAddresses(TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const{

  tree->SetBranchAddress("binNumber", binNumber);
  tree->SetBranchAddress("linkedBins", linkedBins);
  for (int i = 0; i < _dimension; i++) {
    TString lowCornerName  = "lowCorner_"; lowCornerName += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;
    tree->SetBranchAddress(lowCornerName, lowCorner + i);
    tree->SetBranchAddress(highCornerName, highCorner + i);
  }

}

///Load HyperVolumeBinning from a file
///
void HyperVolumeBinning::load(TString filename){

  TFile* file = new TFile(filename, "READ");

  if (file == 0){
    ERROR_LOG << "Could not open TFile in HyperVolumeBinning::load(" << filename << ")";
    return;
  }
  
  loadPrimaryVolumeNumbers(file);

  TTree* tree = (TTree*)file->Get("HyperVolumeBinning");

  if (tree == 0){
    ERROR_LOG << "Could not open TTree in HyperVolumeBinning::load()";
    return;
  }

  //Create branch addresses and link them to TTree
  int binNumber = -1;
  double* lowCorner  = new double [_dimension];
  double* highCorner = new double [_dimension];  
  std::vector<int>* linkedBins = new std::vector<int>();

  setBranchAddresses(tree, &binNumber, lowCorner, highCorner, &linkedBins);
  

  //Loop over the TTree and fill the HyperVolumeBinning
  int nEntries = tree->GetEntries();

  int currentBinNumber = -1;
  HyperVolume* currentHyperVolume = new HyperVolume(_dimension);
  std::vector<int> currentLinkedVolumes;

  //_verbose = true;

  for(int ent = 0; ent < nEntries; ent++){
    tree->GetEntry(ent);
    
    //If the bin number hasn't changed, need to
    //add HyperCube to previous HyperVolume
    if(ent == 0 || currentBinNumber == binNumber){
      currentBinNumber = binNumber;
      currentLinkedVolumes = *linkedBins;
      HyperPoint lowCornerVect (_dimension);
      HyperPoint highCornerVect(_dimension);
      for (int i = 0; i < _dimension; i++){
        lowCornerVect .at(i) = lowCorner [i];
        highCornerVect.at(i) = highCorner[i];
      }
      VERBOSE_LOG << "Adding cuboid to volume";
      currentHyperVolume->addHyperCuboid(lowCornerVect, highCornerVect);
    }

    //if the bin number has changed, need to add
    //previous HyperVolume to HyperVolumeBinning
    //and start a new HyperVolume
    else {
      VERBOSE_LOG << "Adding volume to binning";
      this->addHyperVolume(*currentHyperVolume);
      _linkedHyperVolumes.back() = currentLinkedVolumes;
      delete currentHyperVolume;
      currentHyperVolume = new HyperVolume(_dimension);
      currentLinkedVolumes = *linkedBins;
      currentBinNumber = binNumber;
      HyperPoint lowCornerVect (_dimension);
      HyperPoint highCornerVect(_dimension);
      for (int i = 0; i < _dimension; i++){
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
HyperVolumeBinning::~HyperVolumeBinning(){
  GOODBYE_LOG << "Goodbye from the HyperVolumeBinning() Constructor";
}




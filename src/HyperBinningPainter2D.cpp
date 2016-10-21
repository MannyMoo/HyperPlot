#include "HyperBinningPainter2D.h"

/** Construct a 2D HyperBinningPainter for a given HyperVolumeBinning and HyperPointSet.
This option will draw the HyperPointSet as dots over the binning scheme
*/
HyperBinningPainter2D::HyperBinningPainter2D(BinningBase* binning, HyperPointSet* hyperPoints) :
  HyperBinningPainter(0),
  _binning(binning),
  _hyperPoints(hyperPoints)
{
  if (getBinning().getDimension() != 2){
    ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon";  
  }
  if (hyperPoints != 0){
    if(hyperPoints->getDimension() !=2){
      ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon"; 
    }
  }
  WELCOME_LOG << "Good day from the HyperBinningPainter2D() Constructor";  
}

/** Construct a 2D HyperBinningPainter for a given HyperBinningHistogram */
HyperBinningPainter2D::HyperBinningPainter2D(HyperHistogram* histogram) :
  HyperBinningPainter(histogram),
  _binning(0),
  _hyperPoints(0)
{
  if (getBinning().getDimension() != 2){
    ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon";  
  }
  WELCOME_LOG << "Good day from the HyperBinningPainter2D() Constructor";  
}

/** get the binning  (works for either constructor) */
const BinningBase& HyperBinningPainter2D::getBinning(){
  
  if (_binning != 0) return *_binning;
  return _histogram->getBinning();

}

/** add the bin edges to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawBinEdges(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinEdge(plotter, i);
}

/** add the bin edges to the Plotter (for one HyperVolume) */
void HyperBinningPainter2D::drawBinEdge(RootPlotter2D* plotter, int bin){
  
  for (int i = 0; i < getBinning().getBinHyperVolume(bin).size(); i++){
    HyperCuboid temp (getBinning().getBinHyperVolume(bin).getHyperCuboid(i));
    drawBinEdge(plotter, &temp);
  }
}

/** add the bin edges to the Plotter (for one HyperCuboid) */
void HyperBinningPainter2D::drawBinEdge(RootPlotter2D* plotter, HyperCuboid* bin){

  TLine* leftLine   = new TLine(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getLowCorner().at(0) , bin->getHighCorner().at(1));
  TLine* rightLine  = new TLine(bin->getHighCorner().at(0), bin->getLowCorner().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1));
  TLine* topLine    = new TLine(bin->getLowCorner().at(0) , bin->getHighCorner().at(1), bin->getHighCorner().at(0) , bin->getHighCorner().at(1));
  TLine* bottomLine = new TLine(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getHighCorner().at(0) , bin->getLowCorner().at(1));
  
  leftLine  ->SetLineWidth(1);
  rightLine ->SetLineWidth(1);
  topLine   ->SetLineWidth(1);
  bottomLine->SetLineWidth(1);

  leftLine  ->SetLineColor(kBlack);
  rightLine ->SetLineColor(kBlack);
  topLine   ->SetLineColor(kBlack);
  bottomLine->SetLineColor(kBlack);

//
  leftLine  ->SetLineStyle(3);
  //rightLine ->SetLineStyle(3);
  //topLine   ->SetLineStyle(3);
  bottomLine->SetLineStyle(3);

  plotter->addObject(leftLine);
  //plotter->addObject(rightLine);
  //plotter->addObject(topLine);
  plotter->addObject(bottomLine);

}

/** Code stolen from the THistPainter - choose the colour to
make the bin based on the colour scale of the TH2D */
int HyperBinningPainter2D::getFillColour(double binContents){
  double min = _histogram->getMin();
  double max = _histogram->getMax();
  if (_density == true){
    min = _histogram->getMinDensity();
    max = _histogram->getMaxDensity();
  }

  int ndivz = gStyle->GetNumberContours();
  double dz = max - min;
  double scale = ndivz/dz;
  int ncolors  = gStyle->GetNumberOfColors();

  int color = int(0.01+(binContents-min)*scale);
  int theColor = int((color+0.99)*double(ncolors)/double(ndivz));
  if (theColor > ncolors-1) theColor = ncolors-1;

  return gStyle->GetColorPalette(theColor);
}

/** Draw bin numbers on the all the bins */
void HyperBinningPainter2D::drawBinNumbers(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinNumbers(plotter, i);
}

/** Draw bin number on a single bin */
void HyperBinningPainter2D::drawBinNumbers(RootPlotter2D* plotter, int bin){
  HyperPoint center = getBinning().getBinHyperVolume(bin).getAverageCenter();
  TString label = ""; label += bin;
  plotter->addText( label ,center.at(0), center.at(1), 2, 2, 0.02, 0);

}

/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBins(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawFilledBin(plotter, i);
}

/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBin(RootPlotter2D* plotter, int bin){
  
  double binContent = _histogram->getBinContent(bin);
  if (_density == true ) binContent = _histogram->getFrequencyDensity(bin);

  for (int i = 0; i < getBinning().getBinHyperVolume(bin).size(); i++){
    HyperCuboid temp (getBinning().getBinHyperVolume(bin).getHyperCuboid(i));
    drawFilledBin(plotter, &temp, binContent);
  }

}

/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, double binContents){

  TBox* box = new TBox(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getHighCorner().at(0) , bin->getHighCorner().at(1));
  box->SetFillColor( getFillColour(binContents) );
  box->SetLineWidth(0.0);
  box->SetLineColor( getFillColour(binContents) );
  plotter->addObject(box);

}


/** If HyperPoints provided, add to TH2D */
void HyperBinningPainter2D::addHyperPoints(TH2D* histogram){
  for(unsigned int i = 0; i < _hyperPoints->size(); i++){
    histogram->Fill(_hyperPoints->at(i).at(0), _hyperPoints->at(i).at(1));    
  }  
}

/** Draw the HyperBinningHistogram  */
void HyperBinningPainter2D::draw(TString path){
    
  double x_min = getBinning().getMin(0);
  double x_max = getBinning().getMax(0);
  double y_min = getBinning().getMin(1);
  double y_max = getBinning().getMax(1);
  
  TString xtitle =  _histogram->getNames().getAxisString(0);
  TString ytitle =  _histogram->getNames().getAxisString(1);

  TH2D* histogram = new TH2D("2DHyperBinningPlot", "2DHyperBinningPlot", 100, x_min, x_max, 100, y_min, y_max);
  histogram->GetXaxis()->SetTitle( xtitle );
  histogram->GetYaxis()->SetTitle( ytitle );

  double min = _histogram->getMin();
  double max = _histogram->getMax();
  if (_density == true){
    min = _histogram->getMinDensity();
    max = _histogram->getMaxDensity();
  }

  if (_histogram   != 0){
    //if you want the scale to be correct on the zaxis, we need to fill the
    //histogram with events between min-max (annoying)
    for(int i = 1; i <= 100; i++) for (int j = 0; j <= 100; ++j) histogram->SetBinContent(i, j, max); 
    histogram->SetBinContent(2, 2, min); 
    histogram->SetMaximum(max);
    histogram->SetMinimum(min); 
  } 

  histogram->SetContour(gStyle->GetNumberContours());
  //std::cout << "Am setting number of contours to " << gStyle->GetNumberContours() << std::endl;
  RootPlotter2D plotter(histogram);

  if (_hyperPoints != 0) addHyperPoints(histogram);
  if (_binning     != 0) drawBinEdges  (&plotter);
  if (_histogram   != 0) {
    TBox* box = new TBox(getBinning().getMin(0) , getBinning().getMin(1), getBinning().getMax(0) , getBinning().getMax(1));
    box->SetFillColor( 0 );
    box->SetLineWidth(0.0);
    plotter.addObject(box)  ;
    drawFilledBins(&plotter);
    drawBinEdges(&plotter)  ;
    //drawBinNumbers(&plotter);
  }

  plotter.plot(path, "COLZ");

  delete histogram;
}

/** Constructor  */
HyperBinningPainter2D::~HyperBinningPainter2D(){
  
}

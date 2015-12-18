#include "LHCbStyle.h"

#include "HyperBinningHistogram.h"
#include "LHCbStyle.h"

#include <iostream>
#include <vector>

void testBinning(){
  
  Plotter::s_imageformat = ".pdf";
  
  //You can set whatever dimension you like here - but the plotting only really works for the 1D and 2D cases
  const int dim = 2;
  
  double min = 0.0;
  double max = 1.0;
  
  gRandom->SetSeed(0);
  
  //create the two datasets that we're going to find the chi2 between
  HyperPointSet points1( dim );
  HyperPointSet points2( dim );
  
  //Input the limits of the histogram
  //alternatively, if the different variables have different limits use 
  // HyperCuboid limits( HyperPoint(low0, low1, low2), HyperPoint(up0, up1, up2) );
  
  HyperCuboid limits( dim, min, max );

  for (int i = 0; i < 5000; i++){
    HyperPoint point( dim );
    
    //Fill dataset1 with a randomly generated point
    for (int i = 0; i < dim; i++) point.at(i) = gRandom->Gaus(0.5, 0.25);
    points1.push_back(point);
    
    //Fill dataset2 with a randomly generated point
    for (int i = 0; i < dim; i++) point.at(i) = gRandom->Gaus(0.5, 0.25);
    points2.push_back(point);
  }
  

  
  //Hypername is used to name the variables in a hyper point
  HyperName name("x", "y");

  //Create a histogram with the SMART algorithm using dataset1. This keeps splitting bins
  //so that each of the resulting bins has ~50% of the events. 
  
  std::vector<int> binningDims;    
  binningDims.push_back(0);
  binningDims.push_back(1);

  HyperBinningHistogram hist1(limits, points1, HyperBinningAlgorithms::SMART_MULTI, 
    AlgOption::MinBinContent      (35.0), 
    AlgOption::MinShadowBinContent(35.0), 
    AlgOption::MinBinWidth        (0.0001),
    //AlgOption::UseShadowData      (points2),
    AlgOption::UseWeights         (false),
    AlgOption::RandomSeed         (1),
    AlgOption::StartDimension     (0),
    AlgOption::BinningDimensions  (binningDims),
    AlgOption::DrawAlgorithm("Algorithm")
  );
  
  std::cout << "The maximum bin content is " << hist1.getMax() << std::endl;

  hist1.setNames(name);
  hist1.draw("Histogram1"); 
  hist1.drawDensity("Histogram1_density");  

  //Make a histogram with the same binning and fill it with dataset2

  HyperBinningHistogram hist2( hist1.getBinning() );
  hist2.setNames(name);
  hist2.fill(points2); 
  hist2.draw("Histogram2"); 
  hist2.drawDensity("Histogram2_density");   
  
  //Draw the binning - this doesn't really work at the moment
  
  hist2.getBinning().drawBinning("Binning");
  
  //Find the chi2 between the datasets 
  
  double chi2 = hist1.chi2(hist2);
  int nBins   = hist1.getNBins();
  std::cout << "Chi2 = " << chi2 << "/" << nBins << std::endl;
  
  //Draw a histogram showing the pulls between the histograms

  hist1.pulls(hist2);
  hist1.draw("pullHist");
}



int main(int argc, char **argv) {
  
  LHCbStyle();
  
  testBinning();  

  INFO_LOG << "What's going on";
  ERROR_LOG << "What's going on" << std::endl;

  ERROR_COUNT

}

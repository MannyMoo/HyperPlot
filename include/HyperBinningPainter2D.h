/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class plotting 2D HyperBinningHistograms 
 *
 **/

 
#ifndef HYPERBINNINGPAINTER2D_HH
#define HYPERBINNINGPAINTER2D_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperBinningHistogram.h"
#include "RootPlotter1D.h"
#include "RootPlotter2D.h"
#include "HyperBinningPainter.h"


// Root includes
#include "TH1D.h"
#include "TH2D.h"

// std includes



class HyperBinningPainter2D : public HyperBinningPainter {

  private:
  
  HyperVolumeBinning* _binning; 
  /**< This gets filled in an alternate constructor (usually the HyperBinningHistogram
  get taken from the inhereted HyperBinningPainter class). This allows a HyperPointSet
  to be plotted on top of the HyperVolumeBinning */ 
  HyperPointSet* _hyperPoints;  
  /**< This gets filled in an alternate constructor (usually the HyperBinningHistogram
  get taken from the inhereted HyperBinningPainter class). This allows a HyperPointSet
  to be plotted on top of the HyperVolumeBinning */ 

  const HyperVolumeBinning& getBinning();

  int getFillColour(double binContents);
  
  void addHyperPoints(TH2D* histogram);

  void drawFilledBins(RootPlotter2D* plotter);
  void drawFilledBin(RootPlotter2D* plotter, int bin);
  void drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, double binContents);

  void drawBinEdges(RootPlotter2D* plotter);
  void drawBinEdge(RootPlotter2D* plotter, int bin);
  void drawBinEdge(RootPlotter2D* plotter, HyperCuboid* bin);

  void drawBinNumbers(RootPlotter2D* plotter);

  void drawBinNumbers(RootPlotter2D* plotter, int bin);

  public:

  HyperBinningPainter2D(HyperVolumeBinning* binning, HyperPointSet* hyperPoints = 0);
  HyperBinningPainter2D(HyperBinningHistogram* histogram);

  virtual void draw(TString path = "");

  ~HyperBinningPainter2D();


};



#endif


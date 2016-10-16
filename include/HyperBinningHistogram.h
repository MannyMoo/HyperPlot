/**
  <B>HyperPlot</B>,
  Author: Sam Harnew, sam.harnew@gmail.com ,
  Date: Dec 2015
 

 Binning Algorithms:

  - HyperBinningAlgorithms::SMART             (see HyperBinningMakerSmart for details on the algorithm          )
  - HyperBinningAlgorithms::MINT              (see HyperBinningMakerMint for details on the algorithm           )
  - HyperBinningAlgorithms::MINT_SMART        (see HyperBinningMakerMintSmart for details on the algorithm      )
  - HyperBinningAlgorithms::MINT_RANDOM       (see HyperBinningMakerMintRandomise for details on the algorithm  )
  - HyperBinningAlgorithms::SMART_RANDOM      (see HyperBinningMakerSmartRandomise for details on the algorithm )
  - HyperBinningAlgorithms::LIKELIHOOD        (see HyperBinningMakerLikelihood for details on the algorithm     )
  - HyperBinningAlgorithms::SMART_LIKELIHOOD  (see HyperBinningMakerSmartLikelihood for details on the algorithm)

Binning Algorithm Options:

  - AlgOption::StartDimension     (int dim                  )
  - AlgOption::BinningDimensions  (std::vector<int> dims    )
  - AlgOption::RandomSeed         (int seed                 )
  - AlgOption::MinBinWidth        (double width             )
  - AlgOption::MinBinWidth        (HyperPoint widths        ) 
  - AlgOption::MinBinContent      (double val               )
  - AlgOption::MinShadowBinContent(double val               )
  - AlgOption::UseWeights         (bool   val = true        )
  - AlgOption::UseShadowData      (const HyperPointSet& data)
  - AlgOption::Empty              (                         )
 
 **/
 
#ifndef HYPERBINNINGHISTOGRAM_HH
#define HYPERBINNINGHISTOGRAM_HH

// HyperPlot includes
#include "MessageService.h"
#include "HistogramBase.h"
#include "HyperFunction.h"
#include "HyperVolumeBinning.h"
#include "HyperName.h"

#include "HyperBinningAlgorithms.h"

// Root includes

// std includes

#include <iostream>
#include <fstream>
#include <iomanip>
  
class HyperBinningHistogram : public HistogramBase, public HyperFunction {

  protected:

  HyperVolumeBinning _binning; /**< The HyperVolumeBinning used for the HyperBinningHistogram */
  
  HyperBinningHistogram();
  
  public:
  
  HyperBinningHistogram(
    const HyperCuboid&   binningRange, 
    const HyperPointSet& points, 

    HyperBinningAlgorithms::Alg alg = HyperBinningAlgorithms::MINT, 

    AlgOption opt0 = AlgOption::Empty(),
    AlgOption opt1 = AlgOption::Empty(),  
    AlgOption opt2 = AlgOption::Empty(),
    AlgOption opt3 = AlgOption::Empty(),
    AlgOption opt4 = AlgOption::Empty(),
    AlgOption opt5 = AlgOption::Empty(),
    AlgOption opt6 = AlgOption::Empty(),
    AlgOption opt7 = AlgOption::Empty(),
    AlgOption opt8 = AlgOption::Empty(),
    AlgOption opt9 = AlgOption::Empty()
  );
 
  HyperBinningHistogram(const HyperVolumeBinning& binning);
  HyperBinningHistogram(TString filename);
  HyperBinningHistogram(std::vector<TString> filename);

  
  void setNames( HyperName names ){ _binning.setNames(names); } /**< Set the HyperName (mainly used for axis labels)*/
  HyperName getNames() const {return _binning.getNames();}      /**< Get the HyperName (mainly used for axis labels)*/

  int  fill(const HyperPoint& coords, double weight);
  int  fill(const HyperPoint& coords);
  void fill(const HyperPointSet& points);

  virtual void merge( const HistogramBase& other );

  void merge( TString filenameother );

  
  //Project the HyperBinningHistograms down into one dimension

  void project(TH1D* histogram, const HyperCuboid& cuboid, double content, int dimension) const;
  void project(TH1D* histogram, const HyperVolume& hyperVolume, double content, int dimension) const;
  TH1D project(int dim = 0, int bins = 100, TString name = "projection") const;
  
  void drawProjection    (TString path, int dim = 0, int bins = 100) const;
  void drawAllProjections(TString path, int bins) const;

  void compareProjection    (TString path, int dim, const HyperBinningHistogram& other, int bins = 100) const;
  void compareAllProjections(TString path, const HyperBinningHistogram& other, int bins = 100) const;

  

  HyperBinningHistogram slice(std::vector<int> sliceDims, std::vector<double> sliceVals) const;
  HyperBinningHistogram slice(int dim, double val) const;
  void draw2DSlice   (TString path, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint) const;
  void draw2DSliceSet(TString path, int sliceDimX, int sliceDimY, int sliceSetDim, int nSlices, const HyperPoint& slicePoint) const;
  void draw2DSliceSet(TString path, int sliceDimX, int sliceDimY, int nSlices, const HyperPoint& slicePoint) const;
  void draw2DSliceSet(TString path, int nSlices, const HyperPoint& slicePoint) const;

  HyperCuboid getLimits() const;


  const HyperVolumeBinning& getBinning() const { return _binning; }  /**< get the HyperVolumeBinning */
  
  virtual double getVal(const HyperPoint& point) const;

  virtual double getBinVolume(int bin) const;

  void save(TString filename);

  void load(TString filename);
  
  void setContentsFromFunc(const HyperFunction& func);
  
  void printFull() const;
  
  void saveToTxtFile(TString filename) const;

  void draw(TString path);
  void drawDensity(TString path);
  


  virtual ~HyperBinningHistogram();

};


  //mutable double _nIntegrationsWtrick;
  //mutable double _nIntegrationsWOtrick;
  
  //void printOptimisationStatistics();
  //interpolation stuff
  //HyperPointSet makePointsAtGaussianExtremes(const HyperPoint& mean, const HyperPoint& sigmas, double sigma) const;
  //double intgrateGaussianOverHyperCuboid(const HyperPoint& mean , const HyperPoint& sigmas, const HyperCuboid& cuboid) const;
  //double intgrateGaussianOverHyperVolume(const HyperPoint& point, const HyperPoint& sigmas, const HyperVolume& volume) const;
  //double intgrateGaussianOverBin   (const HyperPoint& point, const HyperPoint& sigmas, int bin) const;
  //double gaussianKernal                 (const HyperPoint& point, const HyperPoint& sigmas) const;
  //HyperPoint findAdaptiveSigma(const HyperPoint& point, const HyperPoint& sigmas) const;
  //std::vector<int> findNHighestContributingKernalBins(const HyperPoint& point, const HyperPoint& sigmas, int n) const;
  //double adaptiveGaussianKernal(const HyperPoint& point, double smoothing = 1.0) const;
  //void reweightDatasetWithAdaptiveGaussianKernal(HyperPointSet& points, double smoothing = 1.0) const;



#endif


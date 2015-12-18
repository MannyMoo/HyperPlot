/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 **/
 
#ifndef HYPERBINNINGMAKER_HH
#define HYPERBINNINGMAKER_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperCuboid.h"
#include "HyperPointSet.h"
#include "HyperVolumeBinning.h"
#include "RootPlotter1D.h"
#include "RootPlotter2D.h"

// Root includes
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class HyperBinningHistogram;


/**
 * HyperBinningMaker is used to adaptively create a HyperVolumeBinning 
 * from a HyperPointSet. This class just contains a suite of useful tools
 * for splitting bins, and creating the bin hierarchy described in HyperVolumeBinning.
 * Specific types of adaptive binning inheret from this class, and are in HyperBinningMakers.h
 *
 * It is also possible to give it a 'shadow HyperPointSet' - this is useful if one is
 * binning the ratio of two HyperPointSet's and require a minium number of events
 * in both the numerator and denominator.
*/
class HyperBinningMaker {

  private:

  protected:


  static bool s_printBinning; /**< print out the status of the binning algorithm */

  //each bin has:
  //   o  a hyper cuboid defining it
  //   o  a set of data which is contained within it
  //   o  a status

  std::vector<HyperCuboid>        _hyperCuboids;          
  /**< vector of HyperCuboid's that defines the binning hierarchy */
  std::vector< std::vector<int> > _linkedBins;           
  /**< each HyperCuboid in the binning hierarchy has linked bins. 
  If these are empty (no linked bins) then this is a true bin. If it 
  contains links, this HyperVolume is just part of the binning hierarchy
  which is used to speed up the binning of events */

  std::vector<HyperPointSet>      _hyperPointSets;        
  /**< records how many HyperPoints in the inital HyperPointSet (that we're going to adaptively bin) 
  fall into each HyperVolume */

  std::vector<HyperPointSet>      _shadowHyperPointSets;  
  /**< records how many HyperPoints in the inital shadow HyperPointSet (that we're going to adaptively bin) 
  fall into each HyperVolume */
  
  /**  Is the volume split as many times as possible, or can I continue to split it  */
  enum VolumeStatus { DONE, CONTINUE };

  std::vector<int>                _status;                
  /**< the status of each HyperVolume i.e. can we continue splitting it into more bins VolumeStatus::CONTINUE 
  or has it been split as many times as possible VolumeStatus::DONE */

  std::vector<int>                _binningDimensions;     
  /**< what dimensions are we allowed to bin in */
  
  bool _shadowAdded;
  /**< has a shadow HyperPointSet been provided  */

  bool _useEventWeights;
  /**< should event weights be considered when creating the binning  */

  double _minimumBinContent;
  /**< the minimum number of events allowed in a bin - if _useEventWeights is true this 
  corresponds to the sum of weights in each bin  */

  double _shadowMinimumBinContent;
  /**< the minimum number of shadow events allowed in a bin - if _useEventWeights is true this 
  corresponds to the sum of weights in each bin  */

  HyperPoint _minimumEdgeLength;
  /**< the minimum bin width in each dimension */

  TRandom* _random;
  /**< random number generator (some binning adaptive binning schemes may require this) */

  bool _drawAlgorithm;
  /**< if this is true, the binning will be drawn after every interation of the algorithm */

  TString _drawAlgorithmDir;
  /**< directory used to draw the binning */

  int _iterationNum;
  /**< what iteration are we on */

  bool isValidBinningDimension(int dimension);

  
  void addBin(const HyperCuboid& hyperCuboid, const HyperPointSet& hyperPointSet, const HyperPointSet& shadowHyperPointSet, int status);
  HyperPointSet filterHyperPointSet(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid, bool print = false) const;
  HyperCuboid splitBelowPoint(int dim, double splitPoint, const HyperCuboid& original) const;
  HyperCuboid splitAbovePoint(int dim, double splitPoint, const HyperCuboid& original) const;

  double getSumOfWeights(const HyperPointSet& hyperPointSet) const;
  double getWeight(const HyperPoint& hyperPoint) const;              //Will return 1 if _useEventWeights is false. If not, get event weight 0.

  double findSmartSplitPoint(int binNumber, int dimension, double dataFraction) const;
  double findSmartSplitPointInt(int binNumber, int dimension, double dataFraction) const;

  double countEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint) const;
  double countEventsInHyperCuboid(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid) const;
  double countShadowEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint) const;

  //int findFirstSplit() const;
  //int resetSplitOrNots();

  double neg2LLH(int binNumber, int dimension, double splitPoint, bool useConstraints = true);
  double nullNeg2LLH(int binNumber);
  
  public:
    
  HyperBinningMaker(const HyperCuboid& binningRange, const HyperPointSet& data);
  
  //The following functions should be called before any binning 
  //binning algorithms commence
  /* ----------------------------------------------------------------*/
  static void setOutputLevel(bool val){s_printBinning = val;}
  /**< set the verbosity of the output - by default this is on */

  void setBinningDimensions(std::vector<int> dims){_binningDimensions = dims;}
  /**< select which dimensions should be binned */

  void addShadowHyperPointSet(const HyperPointSet& data);

  void setSeed(int seed);
  void useEventWeights(bool val = true){_useEventWeights = val;}
  /**< select if weighted event should be used - by default this is off */

  void setMinimumBinContent(double val);
  void setShadowMinimumBinContent(double val);
  void setMinimumEdgeLength(double val);     
  void setMinimumEdgeLength(HyperPoint val);  

  /*----------------------------------------------------------------*/

  int getNumBins() const;
  int getNumHyperVolumes() const;

      

  TH1D* scanSig(int binNumber, int dimension, int nbins, bool useConstraints = true);
  void  getSplitToMinNeg2LLH(double& split, double& sig, int binNumber, int dimension, bool useConstraints = true);

  void  getDimWithLargestSplitSignificance(int& dim, double& split, int binNumber, bool useConstraints = true);


  int likelihoodSplit(int binNumber);
  int likelihoodSplitAll();
  int smartLikelihoodSplit(int binNumber);
  int smartLikelihoodSplitAll();

  //split the bin in a chosen dimension
  int split(int volumeNumber, int dimension, double splitPoint);

  int splitAll(int dimension, double splitPoint);
  int smartSplitAll(int dimension, double dataFraction);
  int smartSplitAllInt(int dimension, double dataFraction);

  //split the bin in a chosen dimension to give a chosen fraction of events in the resulting bin
  int smartSplit(int binNumber, int dimension, double dataFraction);

  int smartMultiSplit(int binNumber, int dimension, int parts);
  int smartMultiSplit(int binNumber, int dimension);
  int smartMultiSplitAll(int dimension);

  int smartSplitInt(int binNumber, int dimension, double dataFraction);

  int smartSplitAllRandomise(double dataFraction = 0.5);
  int splitAllRandomise(double splitPoint = 0.5);


  void drawCurrentState(TString path) const;
  void finishedIteration();
  
  void drawAfterEachIteration(TString path);

  ///Run the binning algorithm 
  virtual void makeBinning() = 0;


  HyperVolumeBinning getHyperVolumeBinning() const;
  HyperBinningHistogram* getHyperBinningHistogram() const;
  HyperBinningHistogram* getShadowHyperBinningHistogram() const;
  HyperBinningHistogram* getRatioHyperBinningHistogram() const;

  virtual ~HyperBinningMaker();

};



#endif


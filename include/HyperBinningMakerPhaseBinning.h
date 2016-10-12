#ifndef HYPERBINNINGMAKERPHASEBINNING_HH
#define HYPERBINNINGMAKERPHASEBINNING_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperBinningMaker.h"
#include "LoadingBar.h"

// Root includes
#include "TMath.h"

// std includes
#include <complex>

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  
 *  
 *  
 *  
 *
 **/
class HyperBinningMakerPhaseBinning : public HyperBinningMaker{
  
  private:
  
  int    _maximumRandWalks;
  int    _numWalkers;
  double     _walkSizeFrac;
  
  int _numberOfCornerSplits;
  int _numberOfPhaseSplits;

  public:
  
  HyperBinningMakerPhaseBinning(const HyperCuboid& binningRange, HyperFunction* func);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  

  int splitByCoord(int volumeNumber, int dimension, HyperPoint& coord);

  int getBinNumFromFunc(HyperPoint& point);


  void walkOrthogonal(HyperPoint& point, HyperCuboid& walkLimit);
  void walk(HyperPoint& point, HyperCuboid& walkLimit);

  HyperPoint getGrad(HyperPoint& point);
  HyperPoint getGradPos(HyperPoint& point, double funcValAtPoint);

  double getSecondDerivative(HyperPoint& point, HyperPoint& vector, double funcValAtPoint, double& deriv);

  int splitDimFromGrad(int volumeNumber, HyperPoint gradient);


  virtual bool passFunctionCriteria(HyperCuboid& cuboid1, HyperCuboid& cuboid2);
  virtual int functionSplitAll(int dimension);

  virtual int functionSplit(int binNumber   , int& dimension);
  int functionSplitRandom  (int volumeNumber, int dimension);
  

  HyperPointSet getSplitCorners( int volumeNumber );
  HyperPointSet getSplitFaces  ( int volumeNumber );
  HyperPointSet getSplitEdges  ( int volumeNumber );

  HyperPoint orderAndTestSplitPoints(HyperPointSet& points, HyperPoint& point, double valAtPoint, HyperPoint gradient);
  int systematicSplit      (int volumeNumber, int dimension, double valAtCenter, HyperPoint gradient);


  int getBinNumFromFuncVal(double phase);
  
  double getLowBinBoundary(int bin);
  
  double getHighBinBoundary(int bin);
  double closestBinBoundary(double val);


  ~HyperBinningMakerPhaseBinning();
  /**< Destructor  */  
};

#endif


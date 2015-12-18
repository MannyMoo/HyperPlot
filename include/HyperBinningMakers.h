#ifndef HYPERBINNINGMAKERS_HH
#define HYPERBINNINGMAKERS_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperBinningMaker.h"

// Root includes

// std includes


/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first splits all bins in the 'startingDim' so that each bin contains 50%
 *  of the HyperPoint%s. 
 *  It then splits all the resulting bins in the dimension 'startingDim + 1' 
 *  using the same method.
 *  This process iterates until the minimum bin content or minimum bin widths 
 *  have been reached.
 *
 **/
class HyperBinningMakerSmart : public HyperBinningMaker {
  
  private:
  
  int _startingDim;  /**< the dimension to start splitting from  */
  
  public:
  
  HyperBinningMakerSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmart();
  /**< Destructor  */  
  
};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first splits all bins in the 'startingDim' so that each bin is half its size. 
 *  It then splits all the resulting bins in the dimension 'startingDim + 1' 
 *  using the same method.
 *  This process iterates until the minimum bin content or minimum bin widths 
 *  have been reached.
 *
 **/
class HyperBinningMakerMint : public HyperBinningMaker{

  private:

  int _startingDim; /**< the dimension to start splitting from  */

  public:

  HyperBinningMakerMint(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMint();
  /**< Destructor  */  
};


/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  \todo describe how this works... 
 *
 **/
class HyperBinningMakerMultiSmart : public HyperBinningMaker {
  
  private:
  
  int _startingDim;  /**< the dimension to start splitting from  */
  
  public:
  
  HyperBinningMakerMultiSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMultiSmart();
  /**< Destructor  */  
  
};




/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first performs the same algorithm as HyperBinningMakerMint, and then
 *  follows this with HyperBinningMakerSmart
 *
 **/
class HyperBinningMakerMintSmart : public HyperBinningMaker{

  private:

  int _startingDim; /**< the dimension to start splitting from  */

  public:

  HyperBinningMakerMintSmart(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMintSmart();
  /**< Destructor  */  
};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first performs the same algorithm as HyperBinningMakerMint, but when
 *  each bin is split, the splitting dimension is chosen at random
 **/
class HyperBinningMakerMintRandomise : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerMintRandomise(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMintRandomise();
  /**< Destructor  */  
};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first performs the same algorithm as HyperBinningMakerSmart, but when
 *  each bin is split, the splitting dimension is chosen at random
 **/
class HyperBinningMakerSmartRandomise : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerSmartRandomise(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmartRandomise();
  /**< Destructor  */  
};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  \todo rememer how this algorithm works - its something along these
 *  lines... Fits the points in the bin with a flat line, then fits it 
 *  with a step fuction (that has 3 free parameters, height before step
 *  height after step, and step point ). The step point that maximises
 *  the likelihood is chosen as the split point. 
 * 
 *  Additionally, you can compare the significance of splitting in different
 *  dimensions using the likelihood ratio of the flat line and the step function.
 *  choose to split in the dimension that gives the biggest significance.
 **/
class HyperBinningMakerLikelihood : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerLikelihood(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerLikelihood();
  /**< Destructor  */  
};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  \todo rememer how this algorithm works - its something along these
 *  lines... follows the same proceedure as HyperBinningMakerLikelihood
 *  but only the dimension is picked using the likelihood method. The
 *  splitting is done using the SMART method (HyperBinningMakerSmart) .
 *
 **/
class HyperBinningMakerSmartLikelihood : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerSmartLikelihood(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmartLikelihood();
  /**< Destructor  */  
};






#endif


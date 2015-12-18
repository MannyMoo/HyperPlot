/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * 
 * This class allowes binning algorithm options to be passed 
 * to a HyperBinningHistogram. An instance of the class can
 * only be made by one of it's static member functions - 
 * for instance:
 *
 * ~~~ {.cpp}
 * AlgOption opt1 = AlgOption::StartDimension(1);
 * AlgOption opt2 = AlgOption::MinBinContent(0.1);
 * AlgOption opt3 = AlgOption::UseWeights(true);
 * ~~~ 
 *
 * These options are inturpreted by the HyperBinningAlgorithms class,
 * which selects the correct binning algorithm and applies all the
 * options passed to it
 **/

#ifndef HYPERBINNINGALGORITHMS_HH
#define HYPERBINNINGALGORITHMS_HH

// HyperPlot includes
#include "HyperBinningMakers.h"
#include "MessageService.h"

// Root includes

// std includes


class AlgOption{ 

  public:

  ///enum containing a list of avalible options one can pass to a 
  ///binning algorithm 
  enum OptionName{   
    EMTPY,                      /**< Empty option */
    START_DIM,                  /**< The dimension to start splitting bins */
    BINNING_DIMS,               /**< The dimensions that the algorithm is allowed to split */
    RAND_SEED,                  /**< Random seed used by the binning algorithm */
    MIN_BIN_WIDTH,              /**< Minimum bin width */
    MIN_BIN_CONTENT,            /**< Minimum bin content */
    MIN_SHADOW_BIN_CONTENT,     /**< Minimum bin content for shadow events */
    USE_WEIGHTS,                /**< Use weights for calculating the bin contents */
    USE_SHADOW_DATA,            /**< Use a show datatset */
    DRAW_ALGORITHM              /**< Draw the binning at each iteration of the algorithm */
  };

  
  private:
  
  OptionName           _optionName;     /**< the option that this particular intstance of the class represents */
  bool                 _bool;           /**< boolean        option */
  int                  _int;            /**< integer        option */
  double               _double;         /**< double         option */
  std::vector<double>  _doublevector;   /**< vector<double> option */
  std::vector<int>     _intvector;      /**< vector<int>    option */
  HyperPoint           _hyperPoint;     /**< HyperPoint     option */
  const HyperPointSet* _hyperPointSet;  /**< HyperPointSet  option */
  TString              _string;         /**< string         option */

  AlgOption();

  public:
  
  static AlgOption Empty              ();
  static AlgOption StartDimension     (int dim );
  static AlgOption BinningDimensions  (std::vector<int> dims);  
  static AlgOption RandomSeed         (int seed);
  static AlgOption MinBinWidth        (double width);
  static AlgOption MinBinWidth        (HyperPoint widths);
  static AlgOption MinBinContent      (double val );
  static AlgOption MinShadowBinContent(double val );
  static AlgOption UseWeights         (bool   val = true);
  static AlgOption UseShadowData      (const HyperPointSet& data);
  static AlgOption DrawAlgorithm      (TString path);

  bool isEmpty();

  OptionName           getOptionName      ();
  bool                 getBoolOpt         ();
  int                  getIntOpt          ();
  double               getDoubleOpt       ();
  TString              getStringOpt       ();
  std::vector<int>     getIntVectorOpt    ();
  std::vector<double>  getDoubleVectorOpt ();
  HyperPoint           getHyperPointOpt   ();
  const HyperPointSet& getHyperPointSetOpt();

  ~AlgOption();

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * This class is used to pick one of the binning algorithms 
 * (in HyperBinningMakers.h) and then sets any additional 
 * options that are passed as AlgOption's
 **/


class HyperBinningAlgorithms{

public:

  /** enum containing a list of the avalible binning algorithms */
  enum Alg{ SMART, MINT, MINT_SMART, MINT_RANDOM ,SMART_RANDOM, LIKELIHOOD, SMART_LIKELIHOOD, SMART_MULTI };


private:

  std::vector<AlgOption> _algOptions;
  /**< list of options that will be given to the HyperBinningMaker */

  Alg _alg;
  /**< the algorithm (HyperBinningMaker) that will be used in a
  particular instance of the class */


  AlgOption getOpt(AlgOption::OptionName name);
  bool optExist(AlgOption::OptionName name);

public:

  HyperBinningAlgorithms(Alg algorithm);

  void addAlgOption(AlgOption option);

  HyperBinningMaker* getHyperBinningMaker(HyperCuboid binningRange, HyperPointSet points);
  ~HyperBinningAlgorithms();

};


#endif


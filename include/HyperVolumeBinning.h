/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * HyperVolumeBinning is a binning scheme where each bin volume is 
 * defined by a HyperVolume.
 *
 **/

/** \class HyperVolumeBinning

Finding the correct bin number is quite a compuationally
slow process if one has to loop over every bin and
check if a HyperPoint falls within that bin volume. It's not unusual to
have millions of HyperPoints that need to be sorted into tens of thousands of bins. 
This would require billions of calulations.
To speed this process up there is the option to build a hierarchy of bins. A schematic below
shows a 1D example. 

~~~ {.cpp}

       HyperVolume Numbers 

 |-------------0-------------| 

 |------1------|------2------| 

 |--3---|---4--|---5---|--6--|

 |-7-|-8| 

           Bin Numbers

 | 0 | 1|   2  |   3   |  4  |

~~~

Imagine we have a HyperPoint that falls into Bin 0. One would first check if it 
falls into HyperVolume 0 (note the distiction here between Bin/HyperVolume Numbers as
indicated by the figure). First we check if it falls into HyperVolume 0, then HyperVolume 1 or 2, then 
HyperVolume 3 or 4, then HyperVolume 7 or 8. 

In this simple example, it took 7 operations, whereas checking each Bin would have taken 
5. Clearly as the number of bins increases, it becomes computationally much less 
expensive to follow this hierarchy approach.

The 

*/



 
#ifndef HYPERVOLUMEBINNING_HH
#define HYPERVOLUMEBINNING_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperPoint.h"
#include "HyperPointSet.h"
#include "HyperCuboid.h"
#include "HyperVolume.h"
#include "RootPlotter1D.h"
#include "RootPlotter2D.h"
#include "HyperName.h"


// Root includes
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class HyperVolumeBinning {

  private:

  void setBranchAddresses   (TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const;
  void createBranches       (TTree* tree, int* binNumber, double* lowCorner, double* highCorner, std::vector<int>** linkedBins) const;
  void saveHyperVolumeToTree(TTree* tree, double* lowCorner, double* highCorner, const HyperVolume& hyperVolume) const;
  
  int followBinLinks(const HyperPoint& coords, int binNumber) const; 

  int _dimension;	/**< dimensionality of the binning */
  
  mutable bool _changed; 
  /**< 
    keep a record to check if the binning changes -
    if so _binNum, _hyperVolumeNumFromBinNum, _averageBinWidth,
    and _minmax need to be redetermined.
  */

  mutable HyperPoint _averageBinWidth;
  /**< store the a Hyperpoint giving the average bin width in each dimension. This
  is calculated once, and reused. */
  mutable HyperCuboid _minmax; /**< store the HyperCuboid that surrounds the binning. This
  is calculated once, and reused. */

  protected:

  HyperName     _names;
  /**< 
    Possible to assign names to each dimension, which will be used for plotting etc
  */  

  std::vector< HyperVolume      > _hyperVolumes;
  /**< 
    std::vector containing all of the HyperVolumes. Usually, not all of these are bins,
    but part of the bin hierarchy that is discussed in the class description.
  */
  
  std::vector< int > _primaryVolumeNumbers;
  /**<
    Usually all bins will be accessed through one primary volume i.e. volume 0 in
    the below example. If one wants to merge two binning schemes together e.g. 
    Binning1 and Binning2 this involes appending the list of HyperVolumes from
    Binning2 to Binning1. This brings a problem when trying to sort a HyperPoint
    into a specific bin i.e. getBinNumber(HyperPoint) - if it belongs to Binning2, 
    it will first have to check that it doesn't fall into any of the volumes in Binning1. 

    To remedy this problem we introduce the list of primary volumes. If this list exists
    i.e. its size is greater than 0, all the primary volumes will be checked first. The
    existance of this list also implies that ALL volumes are linked to the primary volumes.  


           HyperVolume Numbers 
    
     |-------------0-------------| 
    
     |------1------|------2------| 
    
     |--3---|---4--|---5---|--6--|
    
     |-7-|-8| 

  */

  std::vector< std::vector<int> > _linkedHyperVolumes; 
  /**< 
    Every HyperVolume has a _linkedHyperVolumes - if this is empty, this
    means the HyperVolume is a true Bin. If not it is part of the binning hierarchy.
    In the below, HyperVolume 0 would be linked to HyperVolume [1,2], although 
    HyperVolume 4 would be linked to nothing.

    ~~~ {.cpp}
    
           HyperVolume Numbers 
    
     |-------------0-------------| 
    
     |------1------|------2------| 
    
     |--3---|---4--|---5---|--6--|
    
     |-7-|-8| 
    
               Bin Numbers
    
     | 0 | 1|   2  |   3   |  4  |
    
    ~~~

  */

  mutable std::vector< int >              _binNum; 
  /**< 
    Every HyperVolume has a Bin Number, although these will be -1
    if it's not a true bin, and just part of the bin hierarchy.
    These numbers run from 0 to nBins - 1. In the below example
             
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

  mutable std::vector< int >              _hyperVolumeNumFromBinNum;
  /**< 
    The opposite of the member _binNum. This gives the HyperVolume number
    from the Bin Number. In the below example
             
    ~~~ {.cpp}
                                    0   1   2   3   4  
     _hyperVolumeNumFromBinNum = {  7,  8,  4,  5,  6 }


           HyperVolume Numbers 
    
     |-------------0-------------| 
    
     |------1------|------2------| 
    
     |--3---|---4--|---5---|--6--|
    
     |-7-|-8| 
    
               Bin Numbers
    
     | 0 | 1|   2  |   3   |  4  |
    
    ~~~
  */


  void updateCash() const; 
  void updateBinNumbering() const; 
  void updateAverageBinWidth() const;
  void updateMinMax() const;

  public:
  
  HyperVolumeBinning(int dimension);

  //for loading and saving to TFiles
  void load(TString filename);
  void save(TString filename) const;
  void save() const; 
  
  void savePrimaryVolumeNumbers() const;
  void loadPrimaryVolumeNumbers(TFile* file);


  void setNames( HyperName names ){ _names = names; }  /**< set the a name  for each dimension */
  HyperName getNames() const{return _names;}           /**< get the a names for each dimension */
  
  //Used for getting between the bin number and HyperVolume numbers
  int getHyperVolumeNumber(int binNumber) const;
  int getBinNum(int volumeNumber) const;

  void addPrimaryVolumeNumber(int volumeNumber);

  int getNumHyperVolumes() const;  
  int getNumBins() const;
  
  int getBinNum(const HyperPoint& coords) const;

  const int& getDimension () const {return _dimension;}  /**< get the dimensionality of the binning */

  const HyperVolume& getHyperVolume(int volumeNumber) const{return _hyperVolumes.at(volumeNumber);} /**< get one of the HyperVolumes */
  const HyperVolume& getBinHyperVolume(int binNumber) const;

  bool addHyperVolume(const HyperVolume& hyperVolume);
  bool addHyperVolumeLink(int volumeNum, int linkedVolumeNum);

  void mergeBinnings( const HyperVolumeBinning& other );

  HyperPoint getAverageBinWidth() const;

  void drawBin(int dim, int bin, RootPlotter2D& plotter) const;
  void drawBinning(int dim, TString name, TPad* pad = 0, double upperMargin = 0.2, double lowerMargin = 0.2) const;
  void drawBinning(TString name, TPad* pad = 0) const;

  double getMin(int dimension) const;
  double getMax(int dimension) const;
  
  HyperCuboid getLimits() const;

  std::vector<int> findOrderedBinsOnLine( const HyperPoint& point, int dim ) const;
  bool isLineInVolume( const HyperVolume& volume, const HyperPoint& point, int dim ) const;
  bool isLineInCuboid( const HyperCuboid& cuboid, const HyperPoint& point, int dim ) const;

  std::vector<int> findNNearestBins( const HyperPoint& point, int n ) const;
  
  //void print() const;

  ~HyperVolumeBinning();

};



#endif


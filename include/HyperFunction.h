/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * HyperFunction takes a HyperPoint and returns
 * a double. This can be used to reweight HyperPointSets etc.
 *
 **/

 
#ifndef HYPERFUNCTION_HH
#define HYPERFUNCTION_HH

// HyperPlot includes
#include "MessageService.h"
#include "HyperPointSet.h"

// Root includes

// std includes


class HyperFunction{  

  private:

  public:

  HyperFunction(){;} /**< Constructor */

  virtual double getVal(const HyperPoint& point) const = 0; /**< Virtual function that defines a HyperFunction (Map from HyperPoint -> double) */
  
  void reweightDataset(HyperPointSet& points);

  virtual ~HyperFunction(){;} /**< Destructor */

};


#endif
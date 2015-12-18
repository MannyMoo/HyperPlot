#include "HyperFunction.h"

///Reweight a HyperPointSet by the HyperFunction. If weights already
///exist, the existing weights are mulitplied by the HyperFunction
///evaluation. If not, a the HyperFunction evaluation is added as the
///zeroth weight.
void HyperFunction::reweightDataset(HyperPointSet& points){

  int npoints = points.size();

  for (int i = 0; i < npoints; i++){

    HyperPoint& point = points.at(i);

    double val = this->getVal(point);

    int nW = point.numWeights();

    for (int w = 0; w < nW; w++){
      double oldW = point.getWeight(w);
      double newW = oldW*val;
      point.setWeight(w, newW);
    }

    if (nW == 0) point.addWeight(val);

  }


}

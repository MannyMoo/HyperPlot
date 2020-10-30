#ifndef __FOURPIUTILS_H__
#define __FOURPIUTILS_H__

#include <TLorentzVector.h>
#include <include/HyperPoint.h>
#include <utility>

class HyperHistogram;

// Python bindings don't seem to work with functions in namespaces (though it's OK with classes).
// No idea why.
//namespace FourPiUtils{
  double FourPiUtils__getPhi(TLorentzVector particleA1, TLorentzVector particleA2, TLorentzVector particleB1, 
		TLorentzVector particleB2);
  
  double FourPiUtils__getCosTheta(TLorentzVector particle, const TLorentzVector& parent, TLorentzVector grandparent) ;
  
  std::pair<HyperPoint, int>  FourPiUtils__getHyperPoint(const TLorentzVector& particleA1,
							 const TLorentzVector& particleA2,
							 const TLorentzVector& particleB1,
							 const TLorentzVector& particleB2);

  int FourPiUtils__getBinNumber(const HyperHistogram&, const TLorentzVector& particleA1,
				const TLorentzVector& particleA2,
				const TLorentzVector& particleB1, const TLorentzVector& particleB2,
				bool print = false);
  
  void FourPiUtils__useHyperBinning(const std::string&);
//}

#endif

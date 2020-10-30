#include <HyperHistogram.h>
#include <TVector3.h>
#include <FourPiUtils.h>

using namespace std;

//namespace FourPiUtils{
  double FourPiUtils__getPhi(TLorentzVector particleA1, TLorentzVector particleA2,
			     TLorentzVector particleB1, TLorentzVector particleB2) {
  
    //find four-vector of mother particle and boost to its COM frame  
  
    TLorentzVector mother = particleA1 + particleA2 + particleB1 + particleB2;
  
    TVector3 boosttomother = -(mother.BoostVector());
    particleA1.Boost(boosttomother);
    particleA2.Boost(boosttomother);
    particleB1.Boost(boosttomother);
    particleB2.Boost(boosttomother);
  
    //get the three-vectors of each particle
  
    TVector3 particleA1vect = particleA1.Vect();
    TVector3 particleA2vect = particleA2.Vect();
    TVector3 particleB1vect = particleB1.Vect();
    TVector3 particleB2vect = particleB2.Vect();
  
    //imagine the decay proceeds via D -> AB, A->A1A2, B->B1B2
    //find the three-vectors of A and B
  
    TVector3 Avect = ((particleA1 + particleA2).Vect()).Unit();
    TVector3 Bvect = ((particleB1 + particleB2).Vect()).Unit();
  
    //Find the unit vector that is normal to A1 and A2
    TVector3 normAvect = (particleA1vect.Cross(particleA2vect)).Unit();
    //Find the unit vector that is normal to B1 and B2
    TVector3 normBvect = (particleB1vect.Cross(particleB2vect)).Unit();
  
    //find the cosine and sine of the angle phi (angle between decay planes)
    double cosPhi = normAvect.Dot(normBvect);
    double sinPhi = (normAvect.Cross(normBvect)).Dot(Bvect);
  
    //find phi from sin phi and cos phi
    double phi = atan2(sinPhi,cosPhi);
   
    return phi;
  
  }
  
  double FourPiUtils__getCosTheta(TLorentzVector particle, const TLorentzVector& parent,
				  TLorentzVector grandparent) {
    
    //boost the four-vector of the particle and its grandparent (the d meson) to the parent COM frame
    TVector3 boosttoparent = -(parent.BoostVector());
     
    particle.Boost(boosttoparent);
    grandparent.Boost(boosttoparent);
     
    //get the three-vector of the particle and its grandparent
    TVector3 particle3 = particle.Vect();
    TVector3 grandparent3 = grandparent.Vect();
    
    //evaluate cos theta
    double numerator = particle3.Dot(grandparent3);
    double denominator = (particle3.Mag())*(grandparent3.Mag());

    double costhe = numerator/denominator;
     
    return costhe;
  
  }  
  
  void FourPiUtils__useHyperBinning(const string& filename) {
    //The four-vectors that describe a D -> pi+ pi+ pi- pi- decay
    TLorentzVector pip1(17.461242    ,-126.98674   ,-192.13082   ,269.86027    );
    TLorentzVector pip2(551.44241    ,126.98674    ,192.13082    ,613.68429    ); 
    TLorentzVector pim1(-174.56477   ,363.76633    ,0            ,426.94097    ); 
    TLorentzVector pim2(-394.33889   ,-363.76633   ,0            ,554.35448    ); 

    HyperHistogram histogram(filename, "MEMRES READ");  
  
    FourPiUtils__getBinNumber(histogram, pip1, pip2, pim1, pim2, true);
  }
  
  pair<HyperPoint, int> FourPiUtils__getHyperPoint(const TLorentzVector& pip1, const TLorentzVector& pip2,
						   const TLorentzVector& pim1, const TLorentzVector& pim2){
  
    //Find the D meson four-vector
    TLorentzVector d(pip1+pip2+pim1+pim2);
  
    //Find the variables that are used to describe a phase space point
    double mPlus         = (pip1+pip2).M();
    double mMinus        = (pim1+pim2).M();
    double cosThetaPlus  = FourPiUtils__getCosTheta(pip1, pip1+pip2, d);
    double cosThetaMinus = FourPiUtils__getCosTheta(pim1, pim1+pim2, d);
    double phi           = FourPiUtils__getPhi(pip1, pip2, pim1, pim2);
  
    //Use the pion mass and the Dz mass to find the minimum value of
    //mPlus and mMinus possible
    double pionMass = pip1.M();
    //double dzMass   =    d.M();
    double mMin = 2.0*pionMass;
  
    //Find the mPlusPrime and mMinusPrime variables
    double mPlusPrime  = 0.0;
    double mMinusPrime = 0.0;

    if (mMinus > mPlus){
      mPlusPrime  = mPlus  + (mPlus - mMin); 
      mMinusPrime = mMinus + (mPlus - mMin); 
    }
    else {
      mPlusPrime  = mPlus  + (mMinus - mMin); 
      mMinusPrime = mMinus + (mMinus - mMin); 
    }  
  
    //Transform the point so that it falls in the region where
    // cosThetaPlus > 0, cosThetaMinus > 0, phi > 0. Record and 
    // CP flips so that the bin number can also be flipped.
  
  
    if (cosThetaPlus < 0.0){
      cosThetaPlus = -cosThetaPlus;
      phi = phi - TMath::Pi();
    }
    if (cosThetaMinus < 0.0){
      cosThetaMinus = -cosThetaMinus;
      phi = phi- TMath::Pi();
    }
  
    //map phi back into the range [-pi,+pi]
    while (phi < -TMath::Pi()){
      phi += 2.0*TMath::Pi();
    }
    while (phi > TMath::Pi()){
      phi -= 2.0*TMath::Pi();
    }  
  
    double flip = 1;
    if (phi < 0){
      std::swap(cosThetaPlus, cosThetaMinus); 
      std::swap(mPlusPrime  , mMinusPrime  ); 
      phi = -phi;
      flip = -1;
    }  
  
    //The point { mPlusPrime, mMinusPrime, cosThetaPlus, cosThetaMinus, phi } is
    //now guaranteed to fall within the hyper-binning that is used to describe the
    //phase space binning
    HyperPoint psPoint(mPlusPrime, mMinusPrime, cosThetaPlus, cosThetaMinus, phi);
  
    return pair(psPoint, flip);
  }

  int FourPiUtils__getBinNumber(const HyperHistogram& histogram, const TLorentzVector& pip1,
				const TLorentzVector& pip2,
				const TLorentzVector& pim1, const TLorentzVector& pim2, bool print){
    auto pointflip = FourPiUtils__getHyperPoint(pip1, pip2, pim1, pim2);
    auto point = pointflip.first;
    auto flip = pointflip.second;

    int binNumber = histogram.getVal( point );
  
    //if the input phase space point needed to have the CP operation applied, need
    //to multiply it by -1
  
    binNumber *= flip;
  
    if(print){
      //print the output to the screen
      
      std::cout << "The phase space point given (in four-vectors) has coordinates:" << std::endl;
      std::cout << "mPlusPrime    = " << point.at(0) << std::endl;
      std::cout << "mMinusPrime   = " << point.at(1) << std::endl;
      std::cout << "cosThetaPlus  = " << point.at(2) << std::endl;
      std::cout << "cosThetaMinus = " << point.at(3) << std::endl;
      std::cout << "phi           = " << point.at(4) << std::endl;
      std::cout << "CP flip mult  = " << flip        << std::endl;
      
      std::cout << std::endl << "This phase space point falls into bin number " << binNumber << std::endl;
    }
    return binNumber;
  }  
  
//}

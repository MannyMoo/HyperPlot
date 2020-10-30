#ifndef __CISIRESULTS_H__
#define __CISIRESULTS_H__

#include "TMatrixD.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <iomanip>

class CiSiResults{
  
  std::vector<double > _mean;
  std::vector<double > _err;
  std::vector<double > _errp;
  std::vector<double > _errm;
  std::vector<TString> _name;

  TMatrixD _correlations;
    
  public:
  
  CiSiResults(TString filename){
    
    TFile* file = new TFile(filename, "READ");
    if (file == 0){
      std::cout << "TFile named " << filename << " does not exist or could not be loaded" << std::endl;
    }
    
    TTree* tree = dynamic_cast<TTree*>( file->Get("FitParameterSnapshot") );
    if (tree == 0){
      std::cout << "TTree does not exist" << std::endl;
    }
        
    //load the central values and uncertainties of the parameters    
    TString* parameterName = new TString(0);
    double  mean         = 0.0;
    double  err          = 0.0;
    double  errp         = 0.0;
    double  errm         = 0.0;     
    
    tree->SetBranchAddress( "parameterName", &parameterName );
    tree->SetBranchAddress( "mean"         , &mean          );
    tree->SetBranchAddress( "err"          , &err           );   
    tree->SetBranchAddress( "errp"         , &errp          );
    tree->SetBranchAddress( "errm"         , &errm          );      
      
    int nParameters = tree->GetEntries();

    for (int i = 0; i < nParameters; i++){
      tree->GetEntry(i);
  
      _name  .push_back(*parameterName);
      _mean  .push_back(mean         );
      _err   .push_back(err          );
      _errp  .push_back(errp         );
      _errm  .push_back(errm         );
      
    }    
    
    //load the correlation matrix
    _correlations.ResizeTo(nParameters,nParameters);
    _correlations = *dynamic_cast<TMatrixD*>( file->Get("correlationMatrix") );
    
    
  }
  
  int getNumPars(){
    return _name.size();
  }
  
  bool checkParNum(int parNum){
    if ( parNum >= getNumPars() || parNum < 0) {
      std::cout << "Parameter number provided does not exist" << std::endl;
      return false;
    }
    return true;
  }
  
  TString getName(int parNum){
    if ( checkParNum(parNum) == false ) return "ERROR";
    return _name.at(parNum);
  }
  
  double getMean(int parNum){
    if ( checkParNum(parNum) == false ) return -999;
    return _mean.at(parNum);
  }
  
  double getErr(int parNum){
    if ( checkParNum(parNum) == false ) return -999;
    return _err.at(parNum);
  }
  
  double getCorrelation(int parNumi, int parNumj){
    if ( checkParNum(parNumi) == false ) return -999;
    if ( checkParNum(parNumj) == false ) return -999;
    return _correlations(parNumi,parNumj);
  }
  
  void Print(){
    
    std::cout << "The D->4pi hadronic parameters are measured to be:" << std::endl << std::endl;
    
    for (int i = 0; i < getNumPars(); i++){
      std::cout << std::setw(5 ) << std::left << i
                << std::setw(26) << std::left << getName(i) 
                << " = "
                << std::setw(12) << std::left << getMean(i)
                << " ± "
                << std::setw(12) << std::left << getErr(i)
                << std::endl;
    }
    
    std::cout << std::endl << "with correlations:" << std::endl << std::endl;
    std::cout << std::setprecision(2);
    
    _correlations.Print();
    
  }
  
  ~CiSiResults(){
  }

};

#endif

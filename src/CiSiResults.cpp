#include <CiSiResults.h>

int loadresults(){
  
  TString filename = "5binpairs_EqualDeldel/stat.root";
  
  CiSiResults results(filename);
  results.Print();
  
  return 0;
  
}

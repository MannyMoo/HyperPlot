#include "LHCbStyle.h"

#include "HyperBinningHistogram.h"
#include "LHCbStyle.h"

#include <iostream>
#include <vector>


///This class returns a phase at every point in
///two dimensional space. 
class PhaseMotion : public HyperFunction{

  bool _returnBinNum;
  int _opt;  
  int _dim;
  int _nbins;

  public:
  
  PhaseMotion(int dim = 2, int opt = 0, int nbins = 3) :
    _returnBinNum(false), 
    _opt(opt),
    _dim(dim),
    _nbins(nbins)
  {

    if (_opt == 0){
      INFO_LOG << "You have chosen to use phase motion number 0. This " << std::endl;
      INFO_LOG << "only depends on dimensions 0 and 1 (x and y) and is " << std::endl;
      INFO_LOG << "of the form atan2(y, x)" << std::endl;
    }
    else if (_opt == 1){
      INFO_LOG << "You have chosen to use phase motion number 1. This " << std::endl;
      INFO_LOG << "is of the form e^{i pi x} * e^{i pi y} ...  " << std::endl;
    }
    else if (_opt == 2){
      INFO_LOG << "You have chosen to use phase motion number 2. This " << std::endl;
      INFO_LOG << "is of the form e^{i 2 pi r} where r^2 = x^2 + y^2 + ...  " << std::endl;
    }
    else{
      ERROR_LOG << "There is no function number " << _opt << "!" << std::endl;
    }

  }
  
  void returnBinNumber(){_returnBinNum = true; }
  void returnPhase    (){_returnBinNum = false;}


  virtual double getVal(const HyperPoint& point) const{
    
    double phase = 0.0;

    if (_opt == 0){
      double x = point.at(0);
      double y = point.at(1);
      phase = atan2(y, x);
    }
    if (_opt == 1){
      std::complex<double> cval(1.0, 0.0);
      for (int i = 0; i < _dim; i++){
        std::complex<double> ctemp(cos(point.at(i)*TMath::Pi()), sin(point.at(i)*TMath::Pi()));
        cval *= ctemp;
      }      
      phase = arg(cval);
    }
    if (_opt == 2){
      double val = point.norm() * 2.0 * TMath::Pi();
      std::complex<double> cval(cos(val), sin(val));
      phase = arg(cval);
    }
    
    if (_returnBinNum == false){
      return phase;
    }

    int pm = 1;
    if (phase < 0){
      pm = -1;
      phase *= -1;
    }
    
    phase /= TMath::Pi();
    phase *= double(_nbins);
    phase += 1.0;
  
    int bin = floor(phase);
    return bin*pm;
    
  } 


  ~PhaseMotion(){

  }

};

TString GetFunctionBinningDir( int dim, int functionNum, int nbinpairs ){
  
  TString outputdir = "functionBinningExample/";
  gSystem->Exec("mkdir " + outputdir);
  
  outputdir += "dim";
  outputdir += dim;
  outputdir += "_func";
  outputdir += functionNum;  
  outputdir += "_binpairs";
  outputdir += nbinpairs;
  outputdir += "/";
  gSystem->Exec("mkdir " + outputdir);

  return outputdir;

}

TString GetDataBinningDir(int dim ){
  
  TString outputdir = "dataBinningExample/";
  gSystem->Exec("mkdir " + outputdir);
  
  outputdir += "dim";
  outputdir += dim;
  outputdir += "/";
  gSystem->Exec("mkdir " + outputdir);

  return outputdir;

}

void DrawHistogramSlices(TString filename, TString outdir, int dim, int nSlicesPerSet = 50){

  HyperPoint slicePoint(dim, 0.0);
  
  HyperBinningHistogram hist(filename, dim);
  hist.draw2DSliceSet(outdir, nSlicesPerSet, slicePoint);

}

void DataBinningExample(int dim){
  
  TString outputdir = GetDataBinningDir(dim);
  
  double min = 0.0;  
  double max = 1.0;
  
  gRandom->SetSeed(0);
  
  //create the two datasets that we're going to find the chi2 between
  HyperPointSet points1( dim );
  HyperPointSet points2( dim );
  
  //Input the limits of the histogram
  //alternatively, if the different variables have different limits use 
  // HyperCuboid limits( HyperPoint(low0, low1, low2), HyperPoint(up0, up1, up2) );
  
  HyperCuboid limits( dim, min, max );

  for (int i = 0; i < 5000; i++){
    HyperPoint point( dim );
    
    //Fill dataset1 with a randomly generated point
    for (int i = 0; i < dim; i++) point.at(i) = gRandom->Gaus(0.5, 0.25);
    points1.push_back(point);
    
    //Fill dataset2 with a randomly generated point
    for (int i = 0; i < dim; i++) point.at(i) = gRandom->Gaus(0.5, 0.25);
    points2.push_back(point);
  }
  

  //Hypername is used to name the variables in a hyper point
  HyperName name("x", "y");


  //Chose what dimensions you want to bin in. If this option isn't
  //passed, default is to bin in all of them
  std::vector<int> binningDims;    
  for (int i = 0; i < dim; i++){
    binningDims.push_back(i);
  }

  //Create a histogram with the SMART algorithm using dataset1. This keeps splitting bins
  //so that each of the resulting bins has ~50% of the events. 
  
  HyperBinningHistogram hist1(limits, points1, 

    /*** Name of the binning algorithm you want to use     */
    HyperBinningAlgorithms::SMART_MULTI, 

    /***  The minimum number of events allowed in each bin */
    /***  from the HyperPointSet provided (points1)        */
    AlgOption::MinBinContent      (35.0),    

    /*** The minimum number of events allowed in each bin  */
    /*** from the shadow HyperPointSet provided. Providing */
    /*** a shadow set is optional (see option below)       */
    AlgOption::MinShadowBinContent(35.0),    

    /*** This minimum bin width allowed. Can also pass a   */
    /*** HyperPoint if you would like different min bin    */
    /*** widths for each dimension                         */
    AlgOption::MinBinWidth        (0.0001),

    /*** If you want to use a shadow dataset, pass it here.*/
    /*** This is useful when you want to bin the ratio of  */
    /*** two samples.                                      */
    //AlgOption::UseShadowData      (points2),

    /*** If you want to use the sum of weights rather than */
    /*** the number of events, set this to true.           */    
    AlgOption::UseWeights         (false),

    /*** Some algorithms use a random number generator. Set*/
    /*** the seed here                                     */
    AlgOption::RandomSeed         (1),

    /*** What dimesnion would you like to split first? Only*/
    /*** applies to certain algortihms                     */
    AlgOption::StartDimension     (0),

    /*** What dimesnions would you like to bin in?         */
    AlgOption::BinningDimensions  (binningDims),

    /*** Setting this option will make the agorithm draw   */
    /*** the binning scheme at each iteration              */
    AlgOption::DrawAlgorithm(outputdir + "Algorithm")
    
  );
  
  INFO_LOG << "The maximum bin content is " << hist1.getMax() << std::endl;

  hist1.setNames(name);
  hist1.save(outputdir + "Histogram1.root");   
  hist1.draw(outputdir + "Histogram1"); 
  hist1.drawDensity(outputdir + "Histogram1_density");  

  //Make a histogram with the same binning and fill it with dataset2

  HyperBinningHistogram hist2( hist1.getBinning() );
  hist2.setNames(name);
  hist2.fill(points2); 
  hist2.save(outputdir + "Histogram2.root");     
  hist2.draw(outputdir + "Histogram2"); 
  hist2.drawDensity(outputdir + "Histogram2_density");   
  
  // If the dimension is bigger than two, draw slices of the binning.
  if (dim > 2){ 

    int nSlicesPerSet = 5;

    TString slicedir = outputdir + "slices/";
    gSystem->Exec("mkdir " + slicedir);

    DrawHistogramSlices(outputdir + "Histogram1.root", slicedir + "Histogram1", dim, nSlicesPerSet);
    DrawHistogramSlices(outputdir + "Histogram2.root", slicedir + "Histogram2", dim, nSlicesPerSet);    

  }

  //Draw the binning - this doesn't really work at the moment
  
  hist2.getBinning().drawBinning(outputdir + "Binning");
  
  //Find the chi2 between the datasets 
  
  double chi2 = hist1.chi2(hist2);
  int nBins   = hist1.getNBins();
  INFO_LOG << "Chi2 = " << chi2 << "/" << nBins << std::endl;
  
  //Draw a histogram showing the pulls between the histograms

  hist1.pulls(hist2);
  hist1.draw(outputdir + "pullHist");
}


void FunctionBinningExample(int dim, int functionNum, int nbinpairs){
  
  //Get an ouput directory based on the input variables
  TString outputdir = GetFunctionBinningDir(dim, functionNum, nbinpairs);
  
  //define the limits of the binning i.e. [-1, 1] in each dimension
  HyperCuboid   limits(dim, -1.0, 1.0);

  //These are the points we want to bin, but since we're binning a function,
  //we can leave this empty
  HyperPointSet points(dim);

  //A class that inherets from HyperFunction. This maps a point in n-dim
  //space to a real number. In this case, to a phase between -pi and pi.
  //Can select three different functions using functionNum
  PhaseMotion phaseMotion(dim, functionNum, nbinpairs);

  //We want to return the phase rather than the bin number (based on the phase)
  phaseMotion.returnPhase();
  

  //Create a histogram with the FUNC_PHASE algorithm, based on the HyperFunction phaseMotion.

  HyperBinningHistogram hist(limits, points, HyperBinningAlgorithms::FUNC_PHASE,

    /***  The minimum number of events allowed in each bin */
    /***  from the HyperPointSet provided (points)         */   
    AlgOption::MinBinContent      (0.0 ), 

    /*** This minimum bin width allowed. Can also pass a   */
    /*** HyperPoint if you would like different min bin    */
    /*** widths for each dimension                         */
    AlgOption::MinBinWidth        (0.02),

    /*** The hyper function that is used for this algorithm */
    AlgOption::UseFunction        (&phaseMotion),

    /*** Force all bin edges to lie along a grid */
    AlgOption::SnapToGrid         (true),

    /*** The grid (from SnapToGrid) has a granularity of the */
    /*** min bin width multiplied by the grid multiplier */    
    AlgOption::GridMultiplier     (2),
    
    /*** How many bin apirs do you want to split that phase into */
    AlgOption::NumPhaseBinPairs   (nbinpairs)

  );
  
  //make phaseMotion reutrn bin number instead of phase, and set the 
  //contents of the bins using this

  phaseMotion.returnBinNumber();
  hist.setContentsFromFunc(phaseMotion);
  
  //save the histogram to a root file

  hist.save(outputdir + "PhaseBinning.root");
  
  INFO_LOG << "This histogram has " << hist.getNBins() << " bins." << std::endl;
  
  //if the histo is 2 or less dimesnions, draw it. If it's more, draw slices instead
  
  if (dim <= 2){     
    hist.draw(outputdir + "PhaseBinning");
  }
  else{
    TString slicedir = outputdir + "slices/";
    gSystem->Exec("mkdir " + slicedir);    
    DrawHistogramSlices(outputdir + "PhaseBinning.root", slicedir, dim);
  }

}


void PrintHelp(){

  INFO_LOG << "------------ HELP ------------" << std::endl;


  INFO_LOG <<  std::endl;
  std::cout << "--data-binning" << std::endl << std::endl;
  INFO_LOG << "Generate a random n-dimensional dataset, then use a binning alogrithm to " << std::endl;
  INFO_LOG << "split n-dimensional space into volumes that contain approximately the same " << std::endl;
  INFO_LOG << "number of events " << std::endl;

  INFO_LOG <<  std::endl;
  std::cout << "--func-binning" << std::endl << std::endl;
  INFO_LOG << "Run a binning algorithm that is based on a function - specifically," << std::endl;
  INFO_LOG << "a function that maps n-dimensional space to a phase in the range" << std::endl;
  INFO_LOG << "[-pi,+pi]. n-dimensional space is split up into volumes (HyperVolumes) such that the" << std::endl;
  INFO_LOG << "phase variation in each bin is limited to a specified range." << std::endl;
  INFO_LOG << "The range [-pi,+pi] is split into 2N uniform ranges, where N is set" << std::endl;
  INFO_LOG << "using the --bin-pairs option. The physics case is for binned ci si measurements" << std::endl;
  INFO_LOG << "to optimise sensitivity to gamma and charm mixing." << std::endl;

  INFO_LOG <<  std::endl;
  std::cout << "--verbose" << std::endl << std::endl;
  INFO_LOG << "Turn on verbose message service" << std::endl;

  INFO_LOG <<  std::endl;
  std::cout << "--dim" << std::endl << std::endl;
  INFO_LOG << "Choose what dimensionality you would like to use. Works for --func-binning " << std::endl;
  INFO_LOG << "and --data-binning examples. " << std::endl;

  INFO_LOG <<  std::endl;
  std::cout << "--func-num" << std::endl << std::endl;
  INFO_LOG << "Choose what function you would like to choose for the --func-binning example. " << std::endl;
  INFO_LOG << "There are 3 options, 0, 1 or 2.";

  INFO_LOG << std::endl;
  std::cout << "--bin-pairs" << std::endl << std::endl;
  INFO_LOG << "Choose how many bin pairs you would like for the --func-binning example " << std::endl;

  INFO_LOG << std::endl << std::endl;

}

int main(int argc, char** argv) {
  
  LHCbStyle();
  
  //Set the HyperPlot plotter to output in pdf
  Plotter::s_imageformat = ".pdf";

  bool help                   = 0;
  bool dataBinningExample     = 0;
  bool functionBinningExample = 0;
  bool verbose                = 0;

  int nbinpairs    = 3; 
  int functionNum  = 2; 
  int dim          = 2; 

  for(int i = 1; i<argc; i=i+2){
  
    //Options to do with offline selection
    if       (std::string(argv[i])=="--func-binning"   ) { functionBinningExample =  1  ; i--; }
    else if  (std::string(argv[i])=="--data-binning"   ) { dataBinningExample     =  1  ; i--; }
    else if  (std::string(argv[i])=="--help"           ) { help                   =  1  ; i--; }
    else if  (std::string(argv[i])=="--verbose"        ) { verbose                =  1  ; i--; }
    else if  (std::string(argv[i])=="--bin-pairs"      ) { nbinpairs          =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--func-num"       ) { functionNum        =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--dim"            ) { dim                =  atoi(argv[i+1]); }

    else { 
      std::cout << "Entered invalid argument " << argv[i] << std::endl;
      return 0;
    }
  }
  
  if (help) {
    PrintHelp();
    return 0;
  }
  
  if (verbose){
    MessageSerivce::getMessageService()._outputOptions[MessageSerivce::ErrorType::VERBOSE] = true;
  }

  if (dataBinningExample){
    DataBinningExample( dim );
  }
  
  if (functionBinningExample){
    FunctionBinningExample( dim, functionNum, nbinpairs );
  }

  //This will print out how many errors have occured in HyperPlot. 
  ERROR_COUNT
  
  return 0;

}

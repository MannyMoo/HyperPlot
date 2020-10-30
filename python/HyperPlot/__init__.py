
def load_cpp_lib() :
    '''Load the C++ libary for this package into ROOT.'''
    import ROOT
    ROOT.gSystem.Load('libHyperPlotLib.so')
    ROOT.gSystem.Load('libHyperPlotDict.so')

# Uncomment the below if you want to automatically load your C++ libraries
load_cpp_lib()

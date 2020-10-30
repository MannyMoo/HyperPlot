################################################################################
# Package: HyperPlot
################################################################################
gaudi_subdir(HyperPlot v0r1)

gaudi_depends_on_subdirs(GaudiConfUtils)

find_package(PythonLibs)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ROOT COMPONENTS Core RIO Hist MathCore 
		             Graf Graf3d
                	     Postscript Gpad
			     Tree TreePlayer
			     Matrix Physics)
find_package(TBB)

# override the genreflex call to wrap it in the right environment
# set(ROOT_genreflex_CMD ${env_cmd} --xml ${env_xml} ${ROOT_genreflex_CMD})

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
gaudi_add_dictionary(HyperPlot
                     dict/HyperPlotDict.h dict/HyperPlotDict.xml
                     SPLIT_CLASSDEF)

gaudi_add_library(HyperPlotLib
                  src/*.cpp
#                   HyperPlotDict_classdef.cpp
                  PUBLIC_HEADERS include
                  INCLUDE_DIRS Boost ROOT include
                  LINK_LIBRARIES Boost ROOT TBB)

# gaudi_add_module(HyperPlot
#                  src/*.cpp
#                  LINK_LIBRARIES LinkerEvent MCEvent PhysEvent TrackEvent RelationsLib DaVinciKernelLib DaVinciMCKernelLib
# 		 INCLUDE_DIRS Boost ROOT)

# # ensure that genreflex is not called twice (in parallel builds)
# add_dependencies(HyperPlotLib HyperPlotGen)

gaudi_install_python_modules()

#gaudi_add_test(QMTest QMTEST)

gaudi_add_executable(TestHyperBinning example/src/TestBinning.cpp LINK_LIBRARIES HyperPlotLib)

# Error Estimator source files used by targets of project ErrorEstimation
#

file(GLOB all_sources *.cpp *.h)
file(GLOB ReconstLaplacian EELaplacianRecLS.cpp EELaplacianRecLS.h EELaplacianRecWLS.cpp EELaplacianRecWLS.h ErrorEstimator.cpp ErrorEstimator.h)

set(LCCEstimator_SRC ${all_sources} PARENT_SCOPE)
set(LCCEstimatorLaplacian_SRC ${ReconstLaplacian} PARENT_SCOPE)
   
include_directories(.)
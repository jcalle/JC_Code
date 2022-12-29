# common source files used to adaptive classes 
#

file(GLOB sources *.cpp *.h)

set(LCCAdaptive_SRC ${sources} PARENT_SCOPE)
   
include_directories(.)
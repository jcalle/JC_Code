file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers} ${LCCCommon_SRC} ${LCCAdaptive_SRC} ${LCCEstimatorLaplacian_SRC} ${LCCMaterials_SRC})

add_executable(Poisson_H1_ReconstLaplacian_New ${main_SRCS})

target_link_libraries(Poisson_H1_ReconstLaplacian_New ${PZ_LIB})


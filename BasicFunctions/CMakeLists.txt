add_library(BasicFunctions BasicFunctions.cxx EMFunctions.cxx)
target_link_libraries(BasicFunctions PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

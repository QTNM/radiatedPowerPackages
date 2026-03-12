add_library(FieldClasses FieldClasses.cxx FieldPointNR.cxx EMFunctions.cxx)
target_link_libraries(FieldClasses PUBLIC BasicCore ROOTUtils ${ROOT_LIBRARIES})

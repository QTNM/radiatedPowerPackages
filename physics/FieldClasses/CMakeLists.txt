add_library(FieldClasses FieldClasses.cxx FieldPointNR.cxx)
target_link_libraries(FieldClasses PUBLIC BasicCore BasicFunctions ${ROOT_LIBRARIES})

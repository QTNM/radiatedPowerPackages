add_library(BasicFunctions 
    BasicFunctions.cxx 
    TritiumSpectrum.cxx 
)
target_link_libraries(BasicFunctions PUBLIC 
    BasicCore 
    ${ROOT_LIBRARIES}
)

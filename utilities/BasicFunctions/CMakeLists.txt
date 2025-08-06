add_library(BasicFunctions 
    BasicFunctions.cxx 
    PureFunctions.cxx 
    EMFunctions.cxx 
    TritiumSpectrum.cxx 
    ComplexVector3.cxx
)
target_link_libraries(BasicFunctions PUBLIC 
    BasicCore 
    SignalUtils 
    ROOTUtils 
    ${ROOT_LIBRARIES}
)

add_library(ROOTUtils 
    GraphUtils.cxx 
    HistogramUtils.cxx 
    FFTAnalysis.cxx 
    CorrelationAnalysis.cxx
)

target_link_libraries(ROOTUtils PRIVATE 
    BasicCore 
    SignalUtils 
    ${ROOT_LIBRARIES}
)

target_include_directories(ROOTUtils PUBLIC 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)
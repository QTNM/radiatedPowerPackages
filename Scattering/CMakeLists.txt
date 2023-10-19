add_library(Scattering BaseScatter.cxx ElasticScatter.cxx InelasticScatter.cxx)
target_link_libraries(Scattering PRIVATE ${Boost_MATH_LIBRARY})
target_include_directories(Scattering PRIVATE ${Boost_INCLUDE_DIRS})
add_library(Scattering BaseScatter.cxx ElasticScatter.cxx InelasticScatter.cxx)
target_link_libraries(Scattering PRIVATE ${Boost_MATH_LIBRARY})

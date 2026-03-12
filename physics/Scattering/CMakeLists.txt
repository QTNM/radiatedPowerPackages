add_library(Scattering 
  BaseScatter.cxx 
  ElasticScatter.cxx 
  InelasticScatter.cxx
)
target_link_libraries(Scattering PRIVATE 
  BasicCore 
  BasicFunctions
  ${Boost_MATH_LIBRARY}
)

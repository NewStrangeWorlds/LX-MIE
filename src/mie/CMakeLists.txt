add_library(lxmie_lib STATIC mie.cpp size_distribution.cpp)
target_include_directories(lxmie_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(lxmie_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lxmie_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

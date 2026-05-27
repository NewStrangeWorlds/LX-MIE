add_library(lxmie_lib STATIC mie.cpp)
target_include_directories(lxmie_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(lxmie_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

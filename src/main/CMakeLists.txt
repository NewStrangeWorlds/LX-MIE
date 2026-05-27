add_executable(lx_mie main.cpp read_files.cpp write_output.cpp)
target_link_libraries(lx_mie PRIVATE lxmie_lib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lx_mie PRIVATE OpenMP::OpenMP_CXX)
endif()

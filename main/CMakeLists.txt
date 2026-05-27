add_executable(lx_mie main.cpp read_files.cpp write_output.cpp)
target_link_libraries(lx_mie PRIVATE lxmie)

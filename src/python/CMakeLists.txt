find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(lxmie lxmie_bindings.cpp)
target_link_libraries(lxmie PRIVATE lxmie_lib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lxmie PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS lxmie DESTINATION .)

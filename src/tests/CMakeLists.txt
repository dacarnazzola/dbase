# compile test executables
add_executable(test_kinds test_kinds.f90)
target_link_libraries(test_kinds lib_modules)

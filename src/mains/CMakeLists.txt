# compile main executables
add_executable(main main.f90)
target_link_libraries(main lib_modules)

add_executable(ex_particle_filter ex_particle_filter.f90)
target_link_libraries(ex_particle_filter lib_modules)

add_library(CoFlameCore STATIC
    cklib.f
    derf.f
    math.f
    polfit.f
    secondt.f
    tranlib.f
    xerror.f
)
target_link_libraries(CoFlameCore PUBLIC MPI::MPI_Fortran)

add_executable(CoFlame CoFlame.f90)
target_link_libraries(CoFlame CoFlameCore)

add_executable(tran tranfit.f)
target_link_libraries(tran CoFlameCore)

add_executable(chem ckinterp.f)
target_link_libraries(chem CoFlameCore)

install(TARGETS CoFlame chem tran
        RUNTIME DESTINATION bin         COMPONENT Runtime
        LIBRARY DESTINATION lib         COMPONENT Runtime
        ARCHIVE DESTINATION lib/CoFlame COMPONENT Development)


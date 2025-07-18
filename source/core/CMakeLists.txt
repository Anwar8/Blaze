find_package(Eigen3 REQUIRED)
if(WITH_MPI)
    list(APPEND CMAKE_PREFIX_PATH "/Users/anwar/UOE/dissertation/code/Trilinos/trilinos-install")
    find_package(Trilinos REQUIRED)
endif(WITH_MPI)
if(BUILD_STATIC_LIBS)
    file(GLOB CORE_SRC "*.cpp")
    add_library(BlazeCore STATIC ${CORE_SRC})
    if(WITH_MPI)
        target_compile_definitions(BlazeCore PUBLIC WITH_MPI)
        target_include_directories(BlazeCore PUBLIC ${Trilinos_INCLUDE_DIRS})
        target_link_libraries(BlazeCore PUBLIC ${Trilinos_LIBRARIES})
    endif(WITH_MPI)
    if(OMP)
        find_package(OpenMP REQUIRED)
        target_compile_definitions(BlazeCore PUBLIC OMP)
        if(KOKKOS)
            target_compile_definitions(BlazeCore PUBLIC KOKKOS)
            target_include_directories(BlazeCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
            target_link_libraries(BlazeCore PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX Kokkos::kokkos)
        else(KOKKOS)
            target_include_directories(BlazeCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
            target_link_libraries(BlazeCore PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)
        endif(KOKKOS)
    else(OMP)
        if(KOKKOS)
            target_compile_definitions(BlazeCore PUBLIC KOKKOS)
            target_include_directories(BlazeCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
            target_link_libraries(BlazeCore PUBLIC Eigen3::Eigen Kokkos::kokkos)
        else(KOKKOS)
            target_include_directories(BlazeCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
            target_link_libraries(BlazeCore PUBLIC Eigen3::Eigen)
        endif(KOKKOS)
    endif(OMP)   
else(BUILD_STATIC_LIBS)
    add_library(BlazeCore INTERFACE)
    if(WITH_MPI)
        target_compile_definitions(BlazeCore INTERFACE WITH_MPI)
        target_include_directories(BlazeCore INTERFACE ${Trilinos_INCLUDE_DIRS})
        target_link_libraries(BlazeCore INTERFACE ${Trilinos_LIBRARIES})
    endif(WITH_MPI)
    if(OMP)
    find_package(OpenMP REQUIRED)
    target_compile_definitions(BlazeCore INTERFACE OMP)
    if(KOKKOS)
        target_compile_definitions(BlazeCore INTERFACE KOKKOS)
        target_include_directories(BlazeCore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries(BlazeCore INTERFACE Eigen3::Eigen OpenMP::OpenMP_CXX Kokkos::kokkos)
    else(KOKKOS)
        target_include_directories(BlazeCore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries(BlazeCore INTERFACE Eigen3::Eigen OpenMP::OpenMP_CXX)
    endif(KOKKOS)
else(OMP)
    if(KOKKOS)
        target_compile_definitions(BlazeCore INTERFACE KOKKOS)
        target_include_directories(BlazeCore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries(BlazeCore INTERFACE Eigen3::Eigen Kokkos::kokkos)
    else(KOKKOS)
        target_include_directories(BlazeCore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries(BlazeCore INTERFACE Eigen3::Eigen)
    endif(KOKKOS)
endif(OMP)   
endif(BUILD_STATIC_LIBS)



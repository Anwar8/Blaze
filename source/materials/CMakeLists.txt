if(BUILD_STATIC_LIBS)
    file(GLOB MAT_SRC
        "*.cpp"
        "sections/*.cpp"
    )
    add_library(BlazeMaterials STATIC ${MAT_SRC})
    target_include_directories(BlazeMaterials PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/sections)
    target_link_libraries(BlazeMaterials BlazeCore)
else(BUILD_STATIC_LIBS)
    add_library(BlazeMaterials INTERFACE)
    target_include_directories(BlazeMaterials INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/sections)
    target_link_libraries(BlazeMaterials INTERFACE BlazeCore)
endif(BUILD_STATIC_LIBS)

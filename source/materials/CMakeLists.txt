file(GLOB MAT_SRC
    "*.cpp"
    "sections/*.cpp"
)
add_library(BlazeMaterials STATIC ${MAT_SRC})
target_include_directories(BlazeMaterials PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/sections)
target_link_libraries(BlazeMaterials BlazeCore)
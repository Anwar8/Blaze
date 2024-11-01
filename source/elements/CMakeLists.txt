if(BUILD_STATIC_LIBS)
    file(GLOB ELEMENTS_SRC "*.cpp"
        "interface/*.cpp"
        "orientation/*.cpp")
    add_library(BlazeElements STATIC ${ELEMENTS_SRC})
    target_include_directories(BlazeMaterials PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR} 
        ${CMAKE_CURRENT_SOURCE_DIR}/interface
        ${CMAKE_CURRENT_SOURCE_DIR}/orientation)
    target_link_libraries(BlazeElements BlazeCore BlazeMaterials)
else(BUILD_STATIC_LIBS)
    add_library(BlazeElements INTERFACE)
    target_include_directories(BlazeMaterials INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR} 
        ${CMAKE_CURRENT_SOURCE_DIR}/interface
        ${CMAKE_CURRENT_SOURCE_DIR}/orientation)
    target_link_libraries(BlazeElements INTERFACE BlazeCore BlazeMaterials)
endif(BUILD_STATIC_LIBS)

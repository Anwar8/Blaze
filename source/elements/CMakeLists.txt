file(GLOB ELEMENTS_SRC "*.cpp"
    "interface/*.cpp"
    "orientation/*.cpp")

add_library(BlazeElements STATIC ${ELEMENTS_SRC})
target_include_directories(BlazeMaterials PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/interface
    ${CMAKE_CURRENT_SOURCE_DIR}/orientation)
target_link_libraries(BlazeElements BlazeCore BlazeMaterials)
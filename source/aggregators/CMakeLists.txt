file(GLOB AGGREGATORS_SRC "*.cpp")
add_library(BlazeAggregators STATIC ${AGGREGATORS_SRC})
target_include_directories(BlazeAggregators ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(BlazeAggregators BlazeCore BlazeMaterials BlazeElements)

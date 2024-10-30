file(GLOB SOLUTION_SRC "*.cpp")

add_library(BlazeSolution STATIC ${SOLUTION_SRC})
target_include_directories(BlazeSolution PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(BlazeSolution BlazeCore BlazeMaterials BlazeElements BlazeAggregators BlazeManagers)
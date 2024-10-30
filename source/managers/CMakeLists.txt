file(GLOB MANAGERS_SRC
    "boundary_conditions/*.cpp"
    "history/*.cpp"
    "loads/*.cpp"
)
add_library(BlazeManagers STATIC ${MANAGERS_SRC})
target_include_directories(BlazeManagers PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions 
    ${CMAKE_CURRENT_SOURCE_DIR}/history
    ${CMAKE_CURRENT_SOURCE_DIR}/loads)
target_link_libraries(BlazeManagers BlazeCore BlazeMaterials BlazeElements BlazeAggregators)
if(BUILD_STATIC_LIBS)
    file(GLOB MANAGERS_SRC
        "boundary_conditions/*.cpp"
        "history/*.cpp"
        "loads/*.cpp"
        "performance_analysis/*.cpp"
    )
    add_library(BlazeManagers STATIC ${MANAGERS_SRC})
    target_include_directories(BlazeManagers PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions 
        ${CMAKE_CURRENT_SOURCE_DIR}/history
        ${CMAKE_CURRENT_SOURCE_DIR}/loads
        ${CMAKE_CURRENT_SOURCE_DIR}/performance_analysis)
    target_link_libraries(BlazeManagers BlazeCore BlazeMaterials BlazeElements BlazeAggregators)
else(BUILD_STATIC_LIBS)
    add_library(BlazeManagers INTERFACE)
    target_include_directories(BlazeManagers INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions 
        ${CMAKE_CURRENT_SOURCE_DIR}/history
        ${CMAKE_CURRENT_SOURCE_DIR}/loads
        ${CMAKE_CURRENT_SOURCE_DIR}/performance_analysis)
    target_link_libraries(BlazeManagers INTERFACE BlazeCore BlazeMaterials BlazeElements BlazeAggregators)
endif(BUILD_STATIC_LIBS)


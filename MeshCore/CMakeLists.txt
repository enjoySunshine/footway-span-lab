# RecastMeshCore is a dependency-free implementation of the contour-set to
# polygon-mesh step. Recast and DetourTileCache both compile the same header so
# that the offline bake and a tile-cache rebuild cannot drift apart.
add_library(MeshCore INTERFACE)
add_library(RecastNavigation::MeshCore ALIAS MeshCore)

target_include_directories(MeshCore INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>"
)

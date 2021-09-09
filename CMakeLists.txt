cmake_minimum_required(VERSION 3.10)
project(proj_ab C)

set(CMAKE_C_STANDARD 11)

find_library(GDAL_path gdal)
add_library(helper SHARED src/help.c)
target_compile_definitions(helper PUBLIC helper)

add_executable(extract extract.c)
#add_executable(hierarchy hierarchy.c)
add_executable(fix_rasters fix_r_rasters.c)
add_executable(hierarchy_v2 hierarchy_v2.c)

#target_link_libraries(hierarchy PUBLIC "${GDAL_path}" m helper)
target_link_libraries(hierarchy_v2 PUBLIC "${GDAL_path}" m helper)
target_link_libraries(extract PUBLIC "${GDAL_path}" m)
target_link_libraries(fix_rasters PUBLIC "${GDAL_path}" m helper)

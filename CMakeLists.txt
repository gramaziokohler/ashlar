cmake_minimum_required(VERSION 3.16)
project(ashlar)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

# Project Files
set(SOURCES
        src/main.cpp
        src/stone.cpp
        src/four_color.cpp
        src/variational_shape_approximation.cpp
        )

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)

cmake_minimum_required(VERSION 3.24)

project(PathTempo VERSION 0.1.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(CTest)

option(PATH_TEMPO_USE_BUNDLED_DEPENDENCIES "Use the pinned HiGHS and Ruckig submodules" ON)

function(path_tempo_force_o2_in_directory directory)
    get_property(targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
    foreach(target IN LISTS targets)
        get_target_property(imported ${target} IMPORTED)
        get_target_property(type ${target} TYPE)
        get_target_property(o2_enabled ${target} EFFECTIVE_O2_ENABLED)
        if(NOT imported AND NOT o2_enabled
           AND (type STREQUAL "EXECUTABLE" OR type STREQUAL "STATIC_LIBRARY"
                OR type STREQUAL "SHARED_LIBRARY" OR type STREQUAL "MODULE_LIBRARY"
                OR type STREQUAL "OBJECT_LIBRARY"))
            if(MSVC)
                target_compile_options(${target} PRIVATE
                    "$<$<COMPILE_LANGUAGE:C,CXX>:/O2>")
            else()
                target_compile_options(${target} PRIVATE
                    "$<$<COMPILE_LANGUAGE:C,CXX>:-O2>")
            endif()
            set_property(TARGET ${target} PROPERTY EFFECTIVE_O2_ENABLED TRUE)
        endif()
    endforeach()

    get_property(subdirectories DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(subdirectory IN LISTS subdirectories)
        path_tempo_force_o2_in_directory("${subdirectory}")
    endforeach()
endfunction()

function(path_tempo_resolve_dependencies)
    if (PATH_TEMPO_USE_BUNDLED_DEPENDENCIES)
        # Keep implementation dependencies static and avoid HiGHS enabling
        # toolchain-sensitive cross-language LTO for its Unix shared build.
        set(BUILD_SHARED_LIBS OFF)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
    endif()

    if (NOT TARGET highs::highs)
        if (PATH_TEMPO_USE_BUNDLED_DEPENDENCIES)
            set(BUILD_CXX ON)
            set(BUILD_CXX_EXE OFF)
            set(BUILD_TESTING OFF)
            set(BUILD_EXAMPLES OFF)
            set(ZLIB OFF)
            set(HIGHS_NO_DEFAULT_THREADS ON)
            set(HIPO OFF)
            add_subdirectory(third_party/highs EXCLUDE_FROM_ALL)
        else()
            find_package(highs CONFIG REQUIRED)
        endif()
    endif()

    if (NOT TARGET ruckig::ruckig)
        if (PATH_TEMPO_USE_BUNDLED_DEPENDENCIES)
            set(BUILD_EXAMPLES OFF)
            set(BUILD_PYTHON_MODULE OFF)
            set(BUILD_CLOUD_CLIENT OFF)
            set(BUILD_TESTS OFF)
            set(BUILD_BENCHMARK OFF)
            add_subdirectory(third_party/ruckig EXCLUDE_FROM_ALL)

            if (WIN32)
                target_compile_definitions(ruckig PUBLIC _USE_MATH_DEFINES)
            endif()
        else()
            find_package(ruckig CONFIG REQUIRED)
        endif()
    endif()
endfunction()

path_tempo_resolve_dependencies()

add_library(PathTempo STATIC
    src/LinearOptimization.cpp
    src/Planner.cpp
    src/Version.cpp)
add_library(PathTempo::PathTempo ALIAS PathTempo)

target_compile_features(PathTempo PUBLIC cxx_std_23)
target_include_directories(PathTempo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(PathTempo
    PRIVATE
        highs::highs
        ruckig::ruckig)

if (MSVC)
    target_compile_options(PathTempo PRIVATE /W4)
else()
    target_compile_options(PathTempo PRIVATE -Wall -Wextra)
endif()

if (BUILD_TESTING)
    add_executable(path_tempo_tests test/test.cpp)
    target_link_libraries(path_tempo_tests PRIVATE PathTempo::PathTempo)

    if (MSVC)
        target_compile_options(path_tempo_tests PRIVATE /W4)
    else()
        target_compile_options(path_tempo_tests PRIVATE -Wall -Wextra)
    endif()

    add_test(NAME path_tempo_tests COMMAND path_tempo_tests)
endif()

add_executable(path_tempo_adaptive_pockets_example EXCLUDE_FROM_ALL
    examples/adaptive_pockets.cpp
    examples/geometry_text_loader.cpp)
target_link_libraries(path_tempo_adaptive_pockets_example PRIVATE PathTempo::PathTempo)
target_compile_definitions(path_tempo_adaptive_pockets_example PRIVATE
    PATH_TEMPO_EXAMPLE_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples/data")

if (MSVC)
    target_compile_options(path_tempo_adaptive_pockets_example PRIVATE /W4)
else()
    target_compile_options(path_tempo_adaptive_pockets_example PRIVATE -Wall -Wextra)
endif()

# Bundled HiGHS and Ruckig are private implementation dependencies, but their
# Debug defaults would otherwise be -O0. Apply only the optimization option to
# every local C/C++ target so Debug and Release both execute at -O2.
path_tempo_force_o2_in_directory("${CMAKE_CURRENT_SOURCE_DIR}")

configure_package_config_file(
    cmake/PathTempoConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PathTempoConfig.cmake
    INSTALL_DESTINATION lib/cmake/PathTempo)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/PathTempoConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(TARGETS PathTempo
    EXPORT PathTempoTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include)
install(DIRECTORY include/ DESTINATION include)
install(EXPORT PathTempoTargets
    FILE PathTempoTargets.cmake
    NAMESPACE PathTempo::
    DESTINATION lib/cmake/PathTempo)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/PathTempoConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PathTempoConfigVersion.cmake
    DESTINATION lib/cmake/PathTempo)
install(FILES LICENSE DESTINATION share/licenses/PathTempo)

cmake_minimum_required(VERSION 3.24)

project(PathTempo VERSION 0.1.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(CTest)

option(PATH_TEMPO_USE_BUNDLED_DEPENDENCIES "Use the pinned Ruckig submodule" ON)

function(path_tempo_resolve_dependencies)
    if (PATH_TEMPO_USE_BUNDLED_DEPENDENCIES)
        # Keep the implementation dependency static.
        set(BUILD_SHARED_LIBS OFF)
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

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(src/Version.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/generated/Version.cpp" @ONLY)

add_library(PathTempo STATIC
    src/Planner.cpp
    "${CMAKE_CURRENT_BINARY_DIR}/generated/Version.cpp")
add_library(PathTempo::PathTempo ALIAS PathTempo)

target_compile_features(PathTempo PUBLIC cxx_std_23)
target_include_directories(PathTempo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(PathTempo
    PRIVATE
        ruckig::ruckig)

if (MSVC)
    target_compile_options(PathTempo PRIVATE /W4)
else()
    target_compile_options(PathTempo PRIVATE -Wall -Wextra)
endif()

if (BUILD_TESTING)
    add_executable(path_tempo_tests test/test.cpp)
    target_link_libraries(path_tempo_tests PRIVATE PathTempo::PathTempo)
    target_compile_definitions(path_tempo_tests PRIVATE PATH_TEMPO_EXPECTED_VERSION=\"${PROJECT_VERSION}\")

    add_executable(path_tempo_continuous_path_options_tests
        test/continuous_path_options_test.cpp
        examples/continuous_path_geometry.cpp
        examples/continuous_path_options.cpp)
    target_compile_features(path_tempo_continuous_path_options_tests PRIVATE cxx_std_23)

    add_executable(path_tempo_endpoint_exact_arc_tests
        test/endpoint_exact_arc_test.cpp
        examples/endpoint_exact_arc.cpp)
    target_link_libraries(path_tempo_endpoint_exact_arc_tests PRIVATE PathTempo::PathTempo)

    add_executable(path_tempo_geometry_text_loader_tests
        test/geometry_text_loader_test.cpp
        examples/geometry_text_loader.cpp)
    target_compile_features(path_tempo_geometry_text_loader_tests PRIVATE cxx_std_23)

    if (MSVC)
        target_compile_options(path_tempo_tests PRIVATE /W4)
        target_compile_options(path_tempo_continuous_path_options_tests PRIVATE /W4)
        target_compile_options(path_tempo_endpoint_exact_arc_tests PRIVATE /W4)
        target_compile_options(path_tempo_geometry_text_loader_tests PRIVATE /W4)
    else()
        target_compile_options(path_tempo_tests PRIVATE -Wall -Wextra)
        target_compile_options(path_tempo_continuous_path_options_tests PRIVATE -Wall -Wextra)
        target_compile_options(path_tempo_endpoint_exact_arc_tests PRIVATE -Wall -Wextra)
        target_compile_options(path_tempo_geometry_text_loader_tests PRIVATE -Wall -Wextra)
    endif()

    add_test(NAME path_tempo_tests COMMAND path_tempo_tests)
    add_test(NAME path_tempo_continuous_path_options_tests COMMAND path_tempo_continuous_path_options_tests)
    add_test(NAME path_tempo_endpoint_exact_arc_tests COMMAND path_tempo_endpoint_exact_arc_tests)
    add_test(NAME path_tempo_geometry_text_loader_tests COMMAND path_tempo_geometry_text_loader_tests)
endif()

add_executable(path_tempo_continuous_path_example
    examples/continuous_path.cpp
    examples/continuous_path_geometry.cpp
    examples/continuous_path_options.cpp
    examples/endpoint_exact_arc.cpp
    examples/geometry_text_loader.cpp)
target_link_libraries(path_tempo_continuous_path_example PRIVATE PathTempo::PathTempo)

if (MSVC)
    target_compile_options(path_tempo_continuous_path_example PRIVATE /W4)
else()
    target_compile_options(path_tempo_continuous_path_example PRIVATE -Wall -Wextra)
endif()

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

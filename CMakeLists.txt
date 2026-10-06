cmake_minimum_required(VERSION 3.21)

project(geometra
    VERSION 0.1.0
    DESCRIPTION "Librería moderna de geometría computacional en C++20"
    LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Opciones (los extras solo se activan por defecto cuando Geometra es el
# proyecto raíz; al consumirla con add_subdirectory/FetchContent quedan OFF).
# ---------------------------------------------------------------------------
set(GEOMETRA_IS_TOP_LEVEL ${PROJECT_IS_TOP_LEVEL})

option(GEOMETRA_BUILD_TESTS "Compilar los tests" ${GEOMETRA_IS_TOP_LEVEL})
option(GEOMETRA_BUILD_EXAMPLES "Compilar los ejemplos" ${GEOMETRA_IS_TOP_LEVEL})
option(GEOMETRA_BUILD_BENCHMARKS "Compilar los benchmarks" ${GEOMETRA_IS_TOP_LEVEL})
option(GEOMETRA_INSTALL "Generar reglas de instalación" ${GEOMETRA_IS_TOP_LEVEL})
option(GEOMETRA_WARNINGS_AS_ERRORS "Tratar los warnings como errores" OFF)
option(GEOMETRA_SANITIZERS "Compilar con AddressSanitizer y UBSan (GCC/Clang)" OFF)

if(GEOMETRA_IS_TOP_LEVEL)
    # Build Release por defecto en generadores de una sola configuración.
    if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo de build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
    endif()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    if(GEOMETRA_SANITIZERS)
        add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address,undefined)
    endif()
endif()

# C++20 estricto (sin extensiones GNU): entre otras cosas, GCC desactiva así la
# contracción de operaciones en FMA, que rompería los predicados exactos.
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(cmake/CompilerWarnings.cmake)

# ---------------------------------------------------------------------------
# Librería
# ---------------------------------------------------------------------------
configure_file(include/geometra/version.hpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/geometra/version.hpp @ONLY)

add_library(geometra
    src/predicates.cpp
    src/polygon.cpp
    src/convex_hull.cpp
    src/delaunay.cpp)
add_library(geometra::geometra ALIAS geometra)

target_compile_features(geometra PUBLIC cxx_std_20)
target_include_directories(geometra PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Los predicados exactos exigen aritmética IEEE 754 estricta: cada operación en
# double debe redondearse por separado (sin FMA implícito, sin fast-math).
target_compile_options(geometra PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>
    $<$<CXX_COMPILER_ID:MSVC>:/fp:precise>)

geometra_enable_warnings(geometra)

set_target_properties(geometra PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

# ---------------------------------------------------------------------------
# Tests, ejemplos y benchmarks
# ---------------------------------------------------------------------------
if(GEOMETRA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(GEOMETRA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(GEOMETRA_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# ---------------------------------------------------------------------------
# Instalación y paquete CMake (find_package(geometra))
# ---------------------------------------------------------------------------
if(GEOMETRA_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS geometra
        EXPORT geometraTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(DIRECTORY include/geometra
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/geometra/version.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometra)
    install(FILES README.md NOTICE LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})

    set(GEOMETRA_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/geometra)

    install(EXPORT geometraTargets
        NAMESPACE geometra::
        DESTINATION ${GEOMETRA_CMAKE_DIR})

    configure_package_config_file(cmake/geometraConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/geometraConfig.cmake
        INSTALL_DESTINATION ${GEOMETRA_CMAKE_DIR})
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/geometraConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)

    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/geometraConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/geometraConfigVersion.cmake
        DESTINATION ${GEOMETRA_CMAKE_DIR})
endif()

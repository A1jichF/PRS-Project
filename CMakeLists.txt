cmake_minimum_required(VERSION 3.20)
project(PRSErosion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PRS_BUILD_APPS  "Build the erosion_cli program" ON)
option(PRS_BUILD_TESTS "Build the tests"               ON)

# MinGW (MSYS2): statički linkujemo C++ runtime. Inače .exe zavisi od
# libstdc++-6.dll, a ako je na PATH-u prije MSYS2 neka druga verzija te DLL
# (npr. iz Git for Windows), program se sruši odmah pri pokretanju.
if(MINGW)
    add_link_options(-static)
endif()

# ── Biblioteka simulacije (bez ikakve zavisnosti od grafike) ────────────────
add_library(erosion STATIC
    src/Grid.cpp
    src/TerrainGen.cpp
    src/ImageIO.cpp
    src/Diagnostics.cpp
    src/KernelsSeq.cpp
    src/Scenarios.cpp
    src/Simulation.cpp
)
target_include_directories(erosion PUBLIC include)
target_compile_features(erosion PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(erosion PRIVATE /W4 /permissive-)
else()
    target_compile_options(erosion PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ── Programi ────────────────────────────────────────────────────────────────
if(PRS_BUILD_APPS)
    add_executable(erosion_cli apps/erosion_cli/main.cpp)
    target_link_libraries(erosion_cli PRIVATE erosion)
endif()

# ── Testovi ─────────────────────────────────────────────────────────────────
if(PRS_BUILD_TESTS)
    enable_testing()
    add_executable(erosion_tests tests/test_main.cpp tests/test_step1.cpp tests/test_step2_water.cpp tests/test_step4_erosion.cpp
                                 tests/test_step5_transport.cpp)
    target_link_libraries(erosion_tests PRIVATE erosion)
    add_test(NAME erosion_tests COMMAND erosion_tests)
endif()

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
    src/Bedrock.cpp
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

# ── DLL s C interfejsom (include/prs_erosion.h) ─────────────────────────────
# Ovo je ono što koristi renderer: jedan prs_erosion.dll + jedan header, bez
# naših izvornih fajlova u njihovom repozitoriju. Unutra je ista statička
# biblioteka `erosion`, pa su rezultati bit-identični (test_capi.cpp).
# C++ runtime je statički linkovan u DLL (add_link_options(-static) iznad),
# pa DLL ne zavisi od libstdc++-6.dll i sl.
option(PRS_BUILD_DLL "Build prs_erosion.dll (C API)" ON)
if(PRS_BUILD_DLL)
    add_library(prs_erosion SHARED src/ErosionApi.cpp)
    target_link_libraries(prs_erosion PRIVATE erosion)
    target_include_directories(prs_erosion PUBLIC include)
    target_compile_definitions(prs_erosion PRIVATE PRS_EROSION_BUILDING_DLL)
    # prs_erosion.dll umjesto libprs_erosion.dll
    set_target_properties(prs_erosion PROPERTIES
        PREFIX ""
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
    if(MINGW)
        # .def fajl s listom izvezenih funkcija: iz njega se MSVC-ov lib.exe
        # može napraviti prs_erosion.lib, ako engine nekad bude gradio MSVC.
        target_link_options(prs_erosion PRIVATE
            "-Wl,--output-def,${CMAKE_CURRENT_BINARY_DIR}/prs_erosion.def")
    endif()
    if(MSVC)
        target_compile_options(prs_erosion PRIVATE /W4 /permissive-)
    else()
        target_compile_options(prs_erosion PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    # `cmake --install build --prefix <folder>` pravi paket za renderer:
    #   include/prs_erosion.h, bin/prs_erosion.dll, lib/(import lib, .def)
    include(GNUInstallDirs)
    install(TARGETS prs_erosion
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(FILES include/prs_erosion.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    if(MINGW)
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/prs_erosion.def"
            DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
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
                                 tests/test_step5_transport.cpp
                                 tests/test_step6_evaporation.cpp
                                 tests/test_golden.cpp
                                 tests/test_bedrock.cpp
                                 tests/test_river.cpp)
    target_link_libraries(erosion_tests PRIVATE erosion)
    if(PRS_BUILD_DLL)
        # C API testovi idu kroz pravi DLL (ne kroz statičku biblioteku).
        target_sources(erosion_tests PRIVATE tests/test_capi.cpp)
        target_link_libraries(erosion_tests PRIVATE prs_erosion)
    endif()
    add_test(NAME erosion_tests COMMAND erosion_tests)
endif()

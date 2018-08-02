#---Find ROOT using `root-config`---
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOTSYS)
string(STRIP ${ROOTSYS} ROOTSYS)
set(ROOTSYS ${ROOTSYS}/share/root/cmake)
list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Minuit Minuit2 TreePlayer RooFit RooFitCore)

#---Find Boost---
find_package(Boost 1.32.0 COMPONENTS system program_options filesystem regex REQUIRED)
include_directories(${Boost_INCLUDE_DIRS}) 

#---Print some diagnostics---
message(STATUS "CMAKE_PREFIX_PATH:         " ${CMAKE_PREFIX_PATH})
message(STATUS "ROOT_LIBRARIES:            " ${ROOT_LIBRARIES})
message(STATUS "BOOST_ROOT:                " ${BOOST_ROOT})
message(STATUS "BOOST_INCLUDEDIR:          " ${BOOST_INCLUDEDIR})
message(STATUS "BOOST_LIBRARYDIR:          " ${BOOST_LIBRARYDIR})
message(STATUS "Boost_NO_SYSTEM_PATHS:     " ${Boost_NO_SYSTEM_PATHS})
message(STATUS "Boost_INCLUDE_DIRS:        " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARIES:           " ${Boost_LIBRARIES})
message(STATUS "Boost_LIBRARY_DIR:         " ${Boost_LIBRARY_DIR})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Add our own include directory---
include_directories(${oxford-lhcb-wmass_SOURCE_DIR}/include)

#---Add ranges::v3 from git submodule---
include_directories(${oxford-lhcb-wmass_SOURCE_DIR}/range-v3/include)

#---Put our binaries in the bin/ directory---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${oxford-lhcb-wmass_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${oxford-lhcb-wmass_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${oxford-lhcb-wmass_SOURCE_DIR}/lib)

#---Setup our executables---
add_executable(test_options.exe test_options.cxx)
target_link_libraries(test_options.exe LINK_PUBLIC  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable(create_templates.exe create_templates.cxx)
target_link_libraries(create_templates.exe LINK_PUBLIC  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

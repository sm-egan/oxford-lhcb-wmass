#---Add our own include directory---
include_directories(${oxford-lhcb-wmass_SOURCE_DIR}/include)

#---Find ROOT using `root-config` if ${ROOTSYS} isn't set---
# (this is needed on my laptop, but not in the LCG environment)
if(NOT DEFINED ENV{ROOTSYS})
  message(STATUS "Getting ROOTSYS from root-config...")
  execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOTSYS)
  string(STRIP ${ROOTSYS} ROOTSYS)
  list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
endif(NOT DEFINED ENV{ROOTSYS})

#---Setup ROOT---
find_package(ROOT REQUIRED COMPONENTS Minuit Minuit2 TreePlayer RooFit RooFitCore TMVA)

#---Find TBB, but don't make it required---
find_package(TBB)

if(TBB_FOUND)
  add_definitions(-DWITH_TBB) # communicate to the code whether to use TBB or not
  include_directories(${TBB_INCLUDE_DIRS})
  add_definitions(${TBB_DEFINITIONS})
endif(TBB_FOUND)

#---Find Boost---
find_package(Boost 1.32.0 COMPONENTS system program_options filesystem regex REQUIRED)
include_directories(${Boost_INCLUDE_DIRS}) 

#---Print some diagnostics---
message(STATUS "CMAKE_PREFIX_PATH:         " ${CMAKE_PREFIX_PATH})
message(STATUS "ROOT_LIBRARIES:            " ${ROOT_LIBRARIES})
message(STATUS "ROOT_USE_FILE:             " ${ROOT_USE_FILE})
message(STATUS "BOOST_ROOT:                " ${BOOST_ROOT})
message(STATUS "BOOST_INCLUDEDIR:          " ${BOOST_INCLUDEDIR})
message(STATUS "BOOST_LIBRARYDIR:          " ${BOOST_LIBRARYDIR})
message(STATUS "Boost_NO_SYSTEM_PATHS:     " ${Boost_NO_SYSTEM_PATHS})
message(STATUS "Boost_INCLUDE_DIRS:        " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARIES:           " ${Boost_LIBRARIES})
message(STATUS "Boost_LIBRARY_DIR:         " ${Boost_LIBRARY_DIR})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Put our binaries in the bin/ directory---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${oxford-lhcb-wmass_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${oxford-lhcb-wmass_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${oxford-lhcb-wmass_SOURCE_DIR}/lib)

#---Setup our executables---

add_executable(Zsample_analysis.exe Zsample_analysis.cxx)
target_link_libraries(Zsample_analysis.exe LINK_PUBLIC
  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

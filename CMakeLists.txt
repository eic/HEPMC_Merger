#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(SignalBackgroundMerger VERSION 1.0)


## This enables us to use CMAKE_INSTALL_FULL_DATADIR in the code
## (should resolve to ${CMAKE_INSTALL_PREFIX}/share)
include(GNUInstallDirs)
add_definitions(-DCMAKE_INSTALL_FULL_DATADIR="${CMAKE_INSTALL_FULL_DATADIR}")
## Also hand over the source directory for people who don't do make install
add_definitions(-DCMAKE_PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

#----------------------------------------------------------------------------
# Find HepMC3
message ( "---- You can specify a HepMC3 installation using \"-DHepMC3=[path/to/lib]\" " )
# message ( "---- HepMC3=" ${HepMC3} )
find_package(HepMC3 3.2.6 REQUIRED HINTS ${HepMC3} ${HepMC3_DIR} )
message ( "-- HepMC3 library found at ${HEPMC3_LIB}" )
# message ( "-- HepMC3 headers found in ${HEPMC3_INCLUDE_DIR}" )


#----------------------------------------------------------------------------
# Find root, for hepmc.root format
find_package(ROOT 6.20 CONFIG REQUIRED)

#----------------------------------------------------------------------------
set(sources ${PROJECT_SOURCE_DIR}/SignalBackgroundMerger.cpp )
# file( headers )

#----------------------------------------------------------------------------
# Add the executable, and link it to external libraries
#
add_executable(SignalBackgroundMerger ${sources})
target_compile_features(SignalBackgroundMerger PUBLIC cxx_std_17)
target_link_libraries(SignalBackgroundMerger ROOT::Tree ROOT::RIO ROOT::Core HepMC3::HepMC3 HepMC3::rootIO )
target_include_directories(SignalBackgroundMerger
  PRIVATE
  ${HEPMC3_INCLUDE_DIR}
)

  
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SignalBackgroundMerger
  EXPORT SignalBackgroundMerger-export
  DESTINATION bin)

#----------------------------------------------------------------------------
## Generate and install cmake configuration file
install(EXPORT SignalBackgroundMerger-export
  FILE
    SignalBackgroundMerger.cmake
  NAMESPACE
    SignalBackgroundMerger::
  DESTINATION
  cmake
  )




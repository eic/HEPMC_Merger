#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

#----------------------------------------------------------------------------
# Determine project version from the latest git tag (e.g. v2.1.0 -> 2.1.0).
# Falls back to the short commit hash if no tag is found, or "unknown" if
# git is unavailable.
find_package(Git QUIET)
if(GIT_FOUND)
  # Try to get the latest tag
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE _GIT_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _GIT_TAG_RESULT
  )
  if(_GIT_TAG_RESULT EQUAL 0 AND _GIT_TAG MATCHES "^v?[0-9]")
    string(REGEX REPLACE "^v([0-9]+)\\.([0-9]+)\\.?([0-9]*).*" "\\1.\\2.\\3" _PROJECT_VER "${_GIT_TAG}")
    # Strip trailing dot if patch component was absent (e.g. "1.0.")
    string(REGEX REPLACE "\\.$" "" _PROJECT_VER "${_PROJECT_VER}")
    set(_PROJECT_VER_DISPLAY "${_PROJECT_VER}")
  else()
    # No tag found — fall back to short commit hash (not a valid semver,
    # so we use 0.0.0 for CMake's VERSION field and store the hash separately)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE _GIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    set(_PROJECT_VER "0.0.0")
    set(_PROJECT_VER_DISPLAY "0.0.0-${_GIT_HASH}")
  endif()
else()
  set(_PROJECT_VER "0.0.0")
  set(_PROJECT_VER_DISPLAY "unknown")
endif()
message(STATUS "SignalBackgroundMerger version: ${_PROJECT_VER_DISPLAY} (from git)")

project(SignalBackgroundMerger VERSION ${_PROJECT_VER} LANGUAGES CXX)
# Make the full display version (including hash suffix when applicable)
# available as a preprocessor definition and CMake variable
set(PROJECT_VERSION_DISPLAY "${_PROJECT_VER_DISPLAY}")


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
target_link_libraries(SignalBackgroundMerger ROOT::Tree ROOT::RIO ROOT::Core ROOT::Hist HepMC3::HepMC3 HepMC3::rootIO )
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




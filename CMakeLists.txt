#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(SignalBackgroundMerger VERSION 1.0)

#----------------------------------------------------------------------------
# Find HepMC3
message ( "---- You can specify a HepMC3 installation using \"-DHepMC3=[path/to/lib]\" " )
# message ( "---- HepMC3=" ${HepMC3} )
find_package(HepMC3 REQUIRED)
message ( "-- HepMC3 library found at ${HEPMC3_LIB}" )
message ( "-- HepMC3 headers found in ${HEPMC3_INCLUDE_DIR}" )


#----------------------------------------------------------------------------
# Find root, for hepmc.root format
#find_package(ROOT 6.20 CONFIG REQUIRED)

#----------------------------------------------------------------------------
file(GLOB sources "${PROJECT_SOURCE_DIR}/src/*.cpp")
# file( headers )

add_library(deduce_reader_wrapper SHARED include/deduce_reader_wrapper/deduce_reader_wrapper.cpp)
#target_link_libraries(deduce_reader_wrapper ${HEPMC3_LIB})

#----------------------------------------------------------------------------
# Add the executable, and link it to external libraries
#
add_executable(SignalBackgroundMerger ${sources} SignalBackgroundMerger.cpp )
target_compile_features(SignalBackgroundMerger PUBLIC cxx_std_17)
target_link_libraries(SignalBackgroundMerger ${HEPMC3_ROOTIO_LIB} ${HEPMC3_LIB} deduce_reader_wrapper)
target_include_directories(SignalBackgroundMerger
  PRIVATE
  "${HEPMC3_INCLUDE_DIR}"
  include/HEPMC_Merger
  include/argparse
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

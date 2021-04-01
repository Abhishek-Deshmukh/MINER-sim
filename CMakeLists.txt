#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(MinerUC CXX) # Always specify the project type.

set(CMAKE_CXX_STANDARD 14) # Works in recent cmakes
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # Requires CMake 2.8.9+
add_compile_options(-fdiagnostics-color=always)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can sWhy all the casting to pointers? The pointers are temporary and cause errors with rvalue checkinget WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED gdml vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#

# Use modern (CMAKE 3.2-3.4 or so build targets instead of archaic system
# This is quickly becoming the preferred way to do things in CMake 3.5
# This replaces, and is far more flexible than:
include(${Geant4_USE_FILE})
#include(${ROOT_USE_FILE})

# Note that build targets hold ALL the information about build, but
# only need to be added once, in target_link_libraries

# Hopefully the findGeant, etc will add a Geant::Geant and Root::Root target,
# replacing this part

separate_arguments(ROOT_DEFINITIONS)

add_library(geant_lib INTERFACE)
target_include_directories(geant_lib INTERFACE ${Geant4_INCLUDE_DIRS})
target_compile_options(geant_lib INTERFACE ${Geant4_DEFINITIONS})
target_compile_definitions(geant_lib INTERFACE G4USE_STD11) # Needed because of definition bug in 10.2
target_link_libraries(geant_lib INTERFACE ${Geant4_LIBRARIES})


add_library(root_lib INTERFACE)
target_include_directories(root_lib INTERFACE ${ROOT_INCLUDE_DIRS})
target_compile_options(root_lib INTERFACE ${ROOT_DEFINITIONS})
target_link_libraries(root_lib INTERFACE ${ROOT_LIBRARIES})

add_library(CRY INTERFACE)
target_include_directories(CRY INTERFACE /home/mercury/Downloads/CRY/cry_v1.7/src)
target_compile_options(CRY INTERFACE "-no-pie")
set(EXTRA_LIBS ${EXTRA_LIBS} -fPIC -L/home/mercury/Downloads/CRY/cry_v1.7/lib -lCRY)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                  )

file(GLOB cry_sources /home/mercury/Downloads/CRY/cry_v1.7/src/*)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(main_lib ${cry_sources} ${sources} ${headers})

target_include_directories(main_lib PUBLIC include)
target_link_libraries(main_lib root_lib geant_lib CRY)


root_generate_dictionary(
  HistManagerDict            # path to dictionary to generate
  HistManager.hh  # list of includes
  LINKDEF       # list of classes to process
  HistManagerLinkDef.hh                    # ROOT linkDef file
)
add_library(HistManager SHARED HistManagerDict.cxx)
target_link_libraries(HistManager PUBLIC main_lib)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(MINERsim MINER.cc)
target_link_libraries(MINERsim PUBLIC main_lib  HistManager)
# Note: linking to a target is special, it pulls in header dirs, etc. for you.
# Public is a little silly here, since it's an executable
# Also the main library is inherited from HistManger, hit, and track,
# but being clear and reduntant doesn't hurt, since it really does need it
# even without the libs
# Of course, you could also argue that I should include root_lib and geant_lib everywhere...

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Alpsim. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MINERUC_SCRIPTS
    Na22.in
#    MCNP_gamma_pdf.gdat
    vis.mac
    run.mac
    read_gdml.mac
    MINERSetupAugust30.gdml
  )

foreach(_script ${MINERUC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(MINERUC DEPENDS MINERsim)



#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS MINERsim DESTINATION bin)

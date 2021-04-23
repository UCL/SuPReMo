#  ====================================================================================================   
#                                                                                                         
#    SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
#             An implementation of the generalised motion modelling and image registration framework      
#                                                                                                         
#    Copyright (c) University College London (UCL). All rights reserved.                                  
#                                                                                                         
#    This software is distributed WITHOUT ANY WARRANTY; without even                                      
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
#    PURPOSE.                                                                                             
#                                                                                                         
#    See LICENSE.txt in the top level directory for details.                                              
#                                                                                                         
#  ====================================================================================================   




######################################################################
# Set the minimum CMake version.
######################################################################
set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.3)

##################################################################################
# Set some CMake Policies.
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
##################################################################################

set(project_policies
  CMP0057 # Support new ``if()`` IN_LIST operator.
  CMP0063 # Honor visibility properties for all target types.
)
foreach(policy ${project_policies})
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

######################################################################
# Define the project
######################################################################
project(Supremo)

set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake
    ${CMAKE_MODULE_PATH} )

######################################################################
# Define the build add_compile_options
######################################################################
OPTION(USE_OPENMP "To use openMP for multi-CPU processing" ON)
OPTION(BUILD_TESTING "To build the unit tests" OFF)


######################################################################
# Include core dependencies
######################################################################
# OpenMP include
if(USE_OPENMP)
  find_package(OpenMP)
  if(NOT OPENMP_FOUND)
    set(USE_OPENMP OFF CACHE BOOL "To use openMP for multi-CPU processing" FORCE)
    message(WARNING "OpenMP does not appear to be supported by your compiler, forcing USE_OPENMP to OFF")
  else(NOT OPENMP_FOUND)
     message(STATUS "Found OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(NOT OPENMP_FOUND)
endif(USE_OPENMP)

# find the nifti-reg installed libraries
find_package(NiftyReg REQUIRED)
if(NiftyReg_FOUND)
  include_directories( ${NiftyReg_INCLUDE_DIR} )
  message("Found NiftyReg")
else()
  message("Didn't find NiftyReg")
endif()

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Build the core libraries
add_subdirectory(motion-lib)
# Build the extentions if required
add_subdirectory(motion-ext)
# Build the executables
add_subdirectory(motion-apps)
# Build testing if required
if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(motion-test)
endif (BUILD_TESTING)


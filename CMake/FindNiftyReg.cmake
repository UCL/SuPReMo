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


#  File adapted from NifTK 

if (NOT NiftyReg_FOUND)
  find_path( NiftyReg_DIR niftyReg )

  # disabled for now: cuda implementation is currently not supported by supremo
  # so the gpu-related files will never be found.
  if(FALSE AND CUDA_FOUND)

    find_path(NiftyReg_INCLUDE_DIR
      NAME _reg_tools_gpu.h
      PATHS ${NiftyReg_DIR}/include
      NO_DEFAULT_PATH
    )

    find_library(NiftyReg_TOOLS_LIBRARY
      NAMES _reg_tools_gpu _reg_tools_gpu${NiftyReg_DEBUG_POSTFIX}
      PATHS ${NiftyReg_DIR}/lib
      NO_DEFAULT_PATH
    )

  else()

    find_path(NiftyReg_INCLUDE_DIR
      NAME _reg_tools.h
      PATHS ${NiftyReg_DIR}/include
      NO_DEFAULT_PATH
    )

    find_library(NiftyReg_TOOLS_LIBRARY
      NAMES _reg_tools _reg_tools${NiftyReg_DEBUG_POSTFIX}
      PATHS ${NiftyReg_DIR}/lib
      NO_DEFAULT_PATH
    )

  endif()

  if(NiftyReg_TOOLS_LIBRARY AND NiftyReg_INCLUDE_DIR)
    set(NiftyReg_FOUND 1)

    foreach(_library 
            _reg_f3d
            _reg_aladin
            _reg_localTrans
            _reg_femTrans
            _reg_resampling
            _reg_blockMatching
            _reg_globalTrans
            _reg_measure
            _reg_ReadWriteImage
            reg_nifti
            reg_png
            _reg_tools
            _reg_maths
            png
            z)


      find_library( _curLib
                    NAME ${_library} 
                    PATHS ${NiftyReg_DIR}/lib
                    NO_DEFAULT_PATH )

      # Append the found lib
      set(NiftyReg_LIBRARIES ${NiftyReg_LIBRARIES} ${_curLib})
      
      # Need to clear the _curLib variable otherwise find_library will not search again
      unset(_curLib CACHE) 

    endforeach(_library)

    list(REMOVE_DUPLICATES NiftyReg_LIBRARIES)

    get_filename_component( NiftyReg_LIBRARY_DIR ${NiftyReg_TOOLS_LIBRARY} DIRECTORY )

    message( "NiftyReg_INCLUDE_DIR: ${NiftyReg_INCLUDE_DIR}" )
    message( "NiftyReg_LIBRARY_DIR: ${NiftyReg_LIBRARY_DIR}" )
    message( "NiftyReg_LIBRARIES:   ${NiftyReg_LIBRARIES}"   )

  else()
    message( FATAL_ERROR "ERROR: NiftyReg not Found" )
  endif()

endif()

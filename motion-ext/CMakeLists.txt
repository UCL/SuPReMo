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




include_directories(${CMAKE_SOURCE_DIR}/motion-lib)

# The motion-compensatedd image reconstruction
add_library( MoCoReconSuperResolutionIBP
               MoCoReconSuperResolutionIBP.h
               MoCoReconSuperResolutionIBP.cpp )
install(TARGETS MoCoReconSuperResolutionIBP RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES MoCoReconSuperResolutionIBP.h DESTINATION include)

# Define the variable to add dependencies to the executable if required
set(SUPREMO_EXT_LIBS ${SUPREMO_EXT_LIBS} MoCoReconSuperResolutionIBP PARENT_SCOPE)

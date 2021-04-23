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




message( "SUPREMO_EXT_LIBS: ${SUPREMO_EXT_LIBS}" )


include_directories( ${CMAKE_SOURCE_DIR}/motion-lib 
                     ${CMAKE_SOURCE_DIR}/motion-ext )

add_executable( runSupremo runSupremo.cpp )

target_link_libraries( runSupremo  
                         Supremo 
                         CommandLineParser 
                         CorrespondenceModel 
                         ImageSimilarity 
                         ObjectiveFunction
            			 MoCoRecon
            			 Optimiser
            			 ImageAcquisition
                         Transformation
                         GitRevision
                         ${SUPREMO_EXT_LIBS}
                         ${NiftyReg_LIBRARIES} )


set_target_properties(GitRevision PROPERTIES LINKER_LANGUAGE CXX)

# Define the binary install target
install( TARGETS runSupremo
         RUNTIME DESTINATION bin COMPONENT Runtime
         LIBRARY DESTINATION lib COMPONENT Runtime
         ARCHIVE DESTINATION lib COMPONENT Runtime )

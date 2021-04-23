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




# Add support to display the current git revision number
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)
add_library( GitRevision
              GitSHA1.h
              ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp )
install(TARGETS GitRevision RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES GitSHA1.h DESTINATION include)

# Build the libraries
# Supremo
add_library( Supremo 
	             SupremoUtils.h
               Supremo.h 
               Supremo.cpp )
install(TARGETS Supremo RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES SupremoUtils.h Supremo.h DESTINATION include)

# Transformation
add_library( Transformation 
               MathAdditions.h
               BSplineTransformation.h 
               BSplineTransformation.cpp
               SlidingTransformation.h 
               SlidingTransformation.cpp
               Transformation.h
               Transformation.cpp )
install(TARGETS Transformation RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES MathAdditions.h BSplineTransformation.h Transformation.h DESTINATION include)

# The command-line parser
add_library( CommandLineParser 
               CommandLineParser.h 
               CommandLineParser.cpp )
install(TARGETS CommandLineParser RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES CommandLineParser.h DESTINATION include)

# The correspondence model
add_library( CorrespondenceModel 
               CorrespondenceModel.h
               CorrespondenceModel.cpp )
install(TARGETS CorrespondenceModel RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES CorrespondenceModel.h DESTINATION include)

# The image similarity
add_library( ImageSimilarity 
               ImageSimilarity.h 
               SSDImageSimilarity.h
               SSDImageSimilarity.cpp )
install(TARGETS ImageSimilarity RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES SSDImageSimilarity.h DESTINATION include)

# The objective function
add_library( ObjectiveFunction 
               ObjectiveFunction.h
               ObjectiveFunction.cpp )
install(TARGETS ObjectiveFunction RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES ObjectiveFunction.h DESTINATION include)

# The motion-compensatedd image reconstruction
add_library( MoCoRecon
               MoCoReconWeightedAverage.h
               MoCoReconWeightedAverage.cpp
               # MoCoReconSuperResolutionIBP.h
               # MoCoReconSuperResolutionIBP.cpp 
               MoCoRecon.h
               MoCoRecon.cpp )
install(TARGETS MoCoRecon RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES MoCoReconWeightedAverage.h MoCoRecon.h DESTINATION include)

# The optimiser
add_library( Optimiser
                Optimiser.h
                Optimiser.cpp
                ConjugateGradientOptimiser.h
                ConjugateGradientOptimiser.cpp )
install(TARGETS Optimiser RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES Optimiser.h ConjugateGradientOptimiser.h DESTINATION include)


# The image acquisition
add_library( ImageAcquisition
       				ImageAcquisition.h 
      				ImageAcquisition.cpp
      				NoImageAcquisition.h 
      				NoImageAcquisition.cpp
              LowResolutionImageAcquisition.h
              LowResolutionImageAcquisition.cpp)
install(TARGETS ImageAcquisition RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install(FILES ImageAcquisition.h NoImageAcquisition.h LowResolutionImageAcquisition.h DESTINATION include)






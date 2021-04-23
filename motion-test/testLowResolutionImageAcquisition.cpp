// ====================================================================================================   
//                                                                                                        
//   SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
//            An implementation of the generalised motion modelling and image registration framework      
//                                                                                                        
//   Copyright (c) University College London (UCL). All rights reserved.                                  
//                                                                                                        
//   This software is distributed WITHOUT ANY WARRANTY; without even                                      
//   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
//   PURPOSE.                                                                                             
//                                                                                                        
//   See LICENSE.txt in the top level directory for details.                                              
//                                                                                                        
// ====================================================================================================   




#include "_reg_ReadWriteImage.h"
#include "Supremo.h"
#include "CommandLineParser.h"
#include "LowResolutionImageAcquisition.h"
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>



// Tolerances allowed 
constexpr auto EPS_SINGLE = 0.00001;




int main( int argc, char **argv )
{
  // Utilise the command line parser
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-dyn"] = { 1, true, "The dynamic image that defines the geometry of the simulated dynamic image." };
  commandLineOptions["-fullImg"] = { 1, true, "The image in full image space." };
  commandLineOptions["-resSimLow"] = { 1, true, "The expected simulated (dynamic) low resolution image." };


  std::cout << "Made it to here!" << std::endl;

   // Parse the command line
   std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>( argc, argv, commandLineOptions );
   std::cout << parser->getCommandLine() << std::endl;

   // Read the dynamic image (in acquisition space)
   const std::string dynamicImageFileName = parser->getCmdOptionAsString( "-dyn" );
   nifti_image* dynamicImage = nifti_image_read( dynamicImageFileName.c_str(), true );

   // Read the image in full image space (high resolution) from which the acquistiion will be simulated
   const std::string fullResolultionImageFileName = parser->getCmdOptionAsString( "-fullImg" );
   nifti_image* fullResolultionImage = nifti_image_read( fullResolultionImageFileName.c_str(), true );

   // Read the expected simulated low-resulution (simulated dynamic) image
   const std::string simLowResolutionImageFileName = parser->getCmdOptionAsString( "-resSimLow" );
   nifti_image* resultSimulatedLowResImage = nifti_image_read( simLowResolutionImageFileName.c_str(), true );


   // Check if any of the images was not sucessfully read
   if (dynamicImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read the dynamic image: %s", parser->getCmdOptionAsString( "-dyn" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }
   if (fullResolultionImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read image in full image space (high-resolution) image: %s", parser->getCmdOptionAsString( "-fullImg" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }
   if (resultSimulatedLowResImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read the expected result image, simulated low resolution %s", parser->getCmdOptionAsString( "-resSimLow" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }

   // Generate the low-resolution image acquisition object
   std::shared_ptr<LowResolutionImageAcquisition> lowResAcquisition = std::make_shared< LowResolutionImageAcquisition >();
   nifti_image* simulatedLowResolutionImage = lowResAcquisition->SimulateImageAcquisition( fullResolultionImage, dynamicImage );
   return checkImageEquality<float>( simulatedLowResolutionImage, resultSimulatedLowResImage, EPS_SINGLE);
}

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
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>



// Test that checks the correctness of the minimum-sized image in full image space
int main( int argc, char **argv )
{
  // Utilise the command line parser
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-fullImg"] = { 1, true, "The high-resultion image in full image space" };
  commandLineOptions["-dyn"] = { 1, true, "The dynamic image in acquisition space" };
  commandLineOptions["-resMinSizeImg"] = { 1, true, "The expected minimum size image in full image space for the given dynamic image" };

   // Parse the command line
   std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>( argc, argv, commandLineOptions );
   std::cout << parser->getCommandLine() << std::endl;

   // Read the image in full image space (high resolution) from which the acquistiion will be simulated
   const std::string fullResolultionImageFileName = parser->getCmdOptionAsString( "-fullImg" );
   nifti_image* fullResolultionImage = nifti_image_read( fullResolultionImageFileName.c_str(), true );
   
   // Read the reference state image
   const std::string dynamicImageFileName = parser->getCmdOptionAsString( "-dyn" );
   nifti_image* dynamicImage = nifti_image_read( dynamicImageFileName.c_str(), true );

   // Read the expected outcome image that defines the minimum size in full image space on the basis of the dynamic image
   const std::string resultMinSizeImageName = parser->getCmdOptionAsString( "-resMinSizeImg" );
   nifti_image* resultMinSizeImage = nifti_image_read( resultMinSizeImageName.c_str(), true );
   
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
   if (resultMinSizeImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read the expected result image, minimum size in full image space %s", parser->getCmdOptionAsString( "-resMinSizeImg" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }

   // Generate the low-resolution image acquisition object
   std::shared_ptr<LowResolutionImageAcquisition> lowResAcquisition = std::make_shared< LowResolutionImageAcquisition >();

  // Request the minimum-size image in full image space
   nifti_image* generatedMinSizeImage = lowResAcquisition->AllocateMinimumSizeImgInFullImgSpace( fullResolultionImage, dynamicImage, 0 );
   
   bool geometryOK = true;
   if (checkImageGeometryEquality( generatedMinSizeImage, resultMinSizeImage ) == EXIT_FAILURE)
   {
     geometryOK = false;
   }
   
   // Check that the image was properly initialised with zeros
   LowResolutionImageAcquisition::PrecisionType* genImgPtr = (LowResolutionImageAcquisition::PrecisionType*) generatedMinSizeImage->data;
   
   bool voxelContentsOK = true;
   LowResolutionImageAcquisition::PrecisionType zero = static_cast<LowResolutionImageAcquisition::PrecisionType>(0.);
   for (int i=0; i < generatedMinSizeImage->nvox; ++i)
   {
     if (genImgPtr[i] != zero) voxelContentsOK = false;
   }

   // Return based on previous tests
   if (geometryOK && voxelContentsOK)
     return EXIT_SUCCESS;
   else
     return EXIT_FAILURE;
}

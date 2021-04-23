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

constexpr auto EPS_SINGLE = 0.00001;

// Test that checks the correctness of the minimum-sized image in full image space
int main( int argc, char **argv )
{
  // Utilise the command line parser
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-fullImg"] = { 1, true, "The high-resultion image in full image space" };
  commandLineOptions["-dyn"] = { 1, true, "The dynamic image in acquisition space" };
  commandLineOptions["-resAdjImg"] = { 1, true, "The expected adjoint image (i.e. spread out to higher resolution)" };
  commandLineOptions["-resAdjWeightImg"] = { 1, true, "The expected adjoint weights image" };
  commandLineOptions["-saveGenAdjImg"] = { 1, false, "Save the generated image after adjoint" };
  commandLineOptions["-saveGenAdjWeightImg"] = { 1, false, "Save the generated weights image after adjoint" };

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
   const std::string resultAdjointImageName = parser->getCmdOptionAsString( "-resAdjImg" );
   nifti_image* resultAdjointImage = nifti_image_read( resultAdjointImageName.c_str(), true );

   // Read the expected outcome image that defines the minimum size in full image space on the basis of the dynamic image
   const std::string resultAdjointWeightImageName = parser->getCmdOptionAsString( "-resAdjWeightImg" );
   nifti_image* resultAdjointWeightImage = nifti_image_read( resultAdjointWeightImageName.c_str(), true );

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
   if (resultAdjointImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read the expected result adjoint image: %s", parser->getCmdOptionAsString( "-resAdjImg" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }
   if (resultAdjointWeightImage == nullptr)
   {
     char msg[200];
     sprintf_s( msg, "Could not read the expected result adjoint weights image: %s", parser->getCmdOptionAsString( "-resAdjWeightImg" ).c_str() );
     supremo_print_error( msg );
     supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
   }

   // Generate the low-resolution image acquisition object
   std::shared_ptr<LowResolutionImageAcquisition> lowResAcquisition = std::make_shared< LowResolutionImageAcquisition >();

   // Calculate the adjoint image and weights
   lowResAcquisition->CalculateAdjoint( fullResolultionImage, dynamicImage );
   nifti_image* generatedAdjointImage = lowResAcquisition->GetImageAfterAdjoint();
   nifti_image* generatedAdjointWeightsImage = lowResAcquisition->GetWeightsImageAfterAdjoint();

   // Save the output if required
   std::string outAdjointFileName = parser->getCmdOptionAsString( "-saveGenAdjImg" );
   std::string outAdjointWeightsFileName = parser->getCmdOptionAsString( "-saveGenAdjWeightImg" );

   if (!outAdjointFileName.empty())
   {
     reg_io_WriteImageFile( generatedAdjointImage, outAdjointFileName.c_str() );
   }
   if (!outAdjointWeightsFileName.empty())
   {
     reg_io_WriteImageFile( generatedAdjointWeightsImage, outAdjointWeightsFileName.c_str() );
   }


   // Check for equality
   bool adjointImageOK = true;
   bool adjointWeightsImageOK = true;
   
   if (checkImageEquality<float>( generatedAdjointImage, resultAdjointImage, EPS_SINGLE ) == EXIT_FAILURE)
   {
     adjointImageOK = false;
   }
   if (checkImageEquality<float>( generatedAdjointWeightsImage, resultAdjointWeightImage, EPS_SINGLE ) == EXIT_FAILURE)
   {
     adjointWeightsImageOK = false;
   }

   // Return based on previous tests
   if (adjointImageOK && adjointWeightsImageOK)
     return EXIT_SUCCESS;
   else
     return EXIT_FAILURE;
}

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
#include "ImagePyramid.h"
#include "checkEquality.h"
#include <memory>
#include <cmath>


#define EPS 0.000001
#define EPS_SINGLE 0.0001





int main(int argc, char **argv)
{
  // Print some usage information
  if (argc != 6)
  {
    fprintf( stderr, "Usage: %s <inputImage> <numberOfLevels> <numberOfLevelsToPerform> <requestedLevel> <expectedOutputImage>\n", argv[0] );
    return EXIT_FAILURE;
  }

  // Get the input parameters
  char *inputImageName=argv[1];
  int numberOfLevels = atoi( argv[2] );
  int numberOfLevelsToPerform = atoi( argv[3] );
  int requestedLevel = atoi( argv[4] );
  char *expectedImageName = argv[5];
  
  // Read the input image
  nifti_image *inputImage = reg_io_ReadImageFile( inputImageName );
  if (inputImage == NULL)
  {
    supremo_print_error( "The input image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the expected image
  nifti_image *expectedImage = reg_io_ReadImageFile( expectedImageName );
  if (expectedImage == NULL)
  {
    supremo_print_error( "The input image could not be read" );
    return EXIT_FAILURE;
  }

  // Generate the image pyramid
  typedef float VoxelType;
  auto imagePyramid = std::make_shared<ImagePyramid<VoxelType>>();
  imagePyramid->GenerateLevels( inputImage,
                                numberOfLevels,
                                numberOfLevelsToPerform,
                                2, nullptr );

  // Check that when one level too low was requested, NULL is returned
  if (imagePyramid->GetLevel( numberOfLevelsToPerform ) != nullptr)
  {
    supremo_print_error( "Image pyramid return of level above numberOfLevelsToPerform must be NULL" );
    return EXIT_FAILURE;
  }

  if (checkImageEquality<VoxelType>( imagePyramid->GetLevel( requestedLevel ), expectedImage, EPS_SINGLE ) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

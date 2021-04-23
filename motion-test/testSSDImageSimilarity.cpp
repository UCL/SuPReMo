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
#include "SSDImageSimilarity.h"
#include <memory>
#include <cmath>
#include <fstream>


#define EPS 0.000001
#define EPS_SINGLE 0.0001




int main(int argc, char **argv)
{
  // Print some usage information
  if (argc != 4)
  {
    fprintf( stderr, "Usage: %s <inputImage1> <inputImage2> <textFileWithExpecetedSimValue> \n", argv[0] );
    return EXIT_FAILURE;
  }

  // Get the input parameters
  char *inputImageNameA =argv[1];
  char *inputImageNameB = argv[2];
  char *ssdMeasureValTextFile = argv[3];
  
  // Read the input image 1
  nifti_image *inputImageA = reg_io_ReadImageFile( inputImageNameA );
  if (inputImageA == NULL)
  {
    supremo_print_error( "The input image 1 could not be read" );
    return EXIT_FAILURE;
  } 
  
  // Read the input image 2
  nifti_image *inputImageB = reg_io_ReadImageFile( inputImageNameB );
  if (inputImageA == NULL)
  {
    supremo_print_error( "The input image 2 could not be read" );
    return EXIT_FAILURE;
  }

  // Read the expected output
  std::ifstream ssdMeasureFile( ssdMeasureValTextFile, std::ifstream::in );
  
  // Check file was opened correctly
  if (!ssdMeasureFile.is_open())
  {
    char msg[200];
    sprintf_s( msg, "Surrogate signal file could not be opened: %s", ssdMeasureValTextFile );
    supremo_print_error( msg );
    return EXIT_FAILURE;
  }
  float expectedMeasure;
  ssdMeasureFile >> expectedMeasure;


  // Generate the image similarity measure object
  std::shared_ptr<SSDImageSimilarity> ssdMeter = std::make_shared<SSDImageSimilarity>();
  double measuredSSDVal=ssdMeter->GetSimilarityMeasureValueForImages( inputImageA, inputImageB );


  // check that when one level to low was requested, null is returned
  if (fabs( measuredSSDVal - expectedMeasure ) > EPS)
  {
    supremo_print_error( "Measured image similarity measure not as expected" );
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

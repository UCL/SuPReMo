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
#include "BSplineTransformation.h"
#include "SlidingTransformation.h"
#include "CommandLineParser.h"
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>


const double EPS = 0.000001;
const double EPS_SINGLE = 0.001;
const unsigned int MAX_NUM_DEVIATIONS_ALLOWED = 3;



int main(int argc, char **argv)
{
  // Set up the command line parser 
  // Specify a map that lists all known parameters
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-cpg1"] = { 1, true, "Control point grid image region 1" };
  commandLineOptions["-cpg2"] = { 1, true, "Control point grid image region 2" };
  commandLineOptions["-distMap"] = { 1, true, "The image with the signed distance transform that splits the source image into the two regions." };
  commandLineOptions["-source"] = { 1, true, "The source image that is being transformed by the transformation." };
  commandLineOptions["-target"] = { 1, true, "The target image into which space the image will be resampled." };
  commandLineOptions["-expectedDVF"] = { 1, false, "The expected DVF" };
  commandLineOptions["-expectedTransSource"] = { 1, false, "The expected transformed source image" };
  commandLineOptions["-expectedGOCTGrad1"] = { 1, false, "The expected gradient of the gap/overlap constraint term for region 1." };
  commandLineOptions["-expectedGOCTGrad2"] = { 1, false, "The expected gradient of the gap/overlap constraint term for region 2." };


  commandLineOptions["-h"] = { 0, false, "Print a usage message." };
  // Parse the command line
  std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>( argc, argv, commandLineOptions );
  std::cout << parser->getCommandLine() << std::endl;

  // Print usage information
  if (parser->cmdOptionExists( "-h" ))
  {
    fprintf( stderr, "Usage: ToDo" );
    return EXIT_FAILURE;
  }
  
  // Read the first cpg image 
  nifti_image* cpg1Image = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-cpg1" ).c_str() );
  if (cpg1Image == NULL)
  {
    supremo_print_error( "The first cpg image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the second cpg image 
  nifti_image* cpg2Image = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-cpg2" ).c_str() );
  if (cpg2Image == NULL)
  {
    supremo_print_error( "The second cpg image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the signed distance map image 
  nifti_image* distMapImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-distMap" ).c_str() );
  if (distMapImage == NULL)
  {
    supremo_print_error( "The signed distance map image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the target image 
  nifti_image* targetImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-target" ).c_str() );
  if (targetImage == NULL)
  {
    supremo_print_error( "The target image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the source image 
  nifti_image* sourceImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-source" ).c_str() );
  if (sourceImage == NULL)
  {
    supremo_print_error( "The source image could not be read" );
    return EXIT_FAILURE;
  }


  // Read the expected DVF image 
  nifti_image* expectedDVFImage = nullptr;
  if (parser->cmdOptionExists( "-expectedDVF" ))
  {
    expectedDVFImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-expectedDVF" ).c_str() );
    if (expectedDVFImage == NULL)
    {
      supremo_print_error( "The expected DVF image could not be read" );
      return EXIT_FAILURE;
    }
  }

  // Read the expected warped source image 
  nifti_image* expectedTransformedSourceImage = nullptr;
  if (parser->cmdOptionExists( "-expectedTransSource" ))
  {
    expectedTransformedSourceImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-expectedTransSource" ).c_str() );
    if (expectedTransformedSourceImage == NULL)
    {
      supremo_print_error( "The expected transformed source image could not be read" );
      return EXIT_FAILURE;
    }
  }


  // Read the expected gradients if required
  nifti_image* expectedGOCTGrad1Image = nullptr;
  if (parser->cmdOptionExists( "-expectedGOCTGrad1" ))
  {
    expectedGOCTGrad1Image = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-expectedGOCTGrad1" ).c_str() );

    if (expectedGOCTGrad1Image == nullptr)
    {
      supremo_print_error( "The expected gap/overlap constraint gradient image could not be read" );
      return EXIT_FAILURE;
    }
  }

  // Read the expected gradients if required
  nifti_image* expectedGOCTGrad2Image = nullptr;
  if (parser->cmdOptionExists( "-expectedGOCTGrad2" ))
  {
    expectedGOCTGrad2Image = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-expectedGOCTGrad2" ).c_str() );

    if (expectedGOCTGrad2Image == nullptr)
    {
      supremo_print_error( "The expected gap/overlap constraint gradient image could not be read" );
      return EXIT_FAILURE;
    }
  }


  // Get the cpg spacing from the image information
  float* cpg1Spacing = new float[cpg1Image->ndim];
  float* cpg2Spacing = new float[cpg2Image->ndim];

  // Dimensions of the different CPGs should be the same
  if (cpg1Image->nu >= 2)
  {
    cpg1Spacing[0] = cpg1Image->dx;
    cpg1Spacing[1] = cpg1Image->dy;
  }
  if (cpg1Image->nu > 2)
  {
    cpg1Spacing[2] = cpg1Image->dz;
  }
  if (cpg2Image->nu >= 2)
  {
    cpg2Spacing[0] = cpg2Image->dx;
    cpg2Spacing[1] = cpg2Image->dy;
  }
  if (cpg2Image->nu > 2)
  {
    cpg2Spacing[2] = cpg2Image->dz;
  }

  // Generate the b-spline transformations
  std::shared_ptr<BSplineTransformation> bsplTrafo1 = std::make_shared<BSplineTransformation>( targetImage, 1, cpg1Spacing );
  std::shared_ptr<BSplineTransformation> bsplTrafo2 = std::make_shared<BSplineTransformation>( targetImage, 1, cpg2Spacing );
   
  // Check that the generated B-spline transformation and the input cpg have the same number of parameters
  if (bsplTrafo1->GetNumberOfParameters() != cpg1Image->nvox)
  {
    supremo_print_error( "Wrong number of parameters in BSpline-Transformation 1" );
    return EXIT_FAILURE;
  }
  if (bsplTrafo2->GetNumberOfParameters() != cpg2Image->nvox)
  {
    supremo_print_error( "Wrong number of parameters in BSpline-Transformation 2" );
    return EXIT_FAILURE;
  }

  // Set the constraints to zero in the first place
  bsplTrafo1->SetBendingEnergyWeight( 0.0 );
  bsplTrafo2->SetBendingEnergyWeight( 0.0 );
  bsplTrafo1->SetLinearEnergyWeight( 0.0 );
  bsplTrafo2->SetLinearEnergyWeight( 0.0 );
  
  // Create the sliding transformation from both bspline transformations
  std::shared_ptr<SlidingTransformation> slidingTrafo = std::make_shared<SlidingTransformation>( bsplTrafo1, bsplTrafo2, targetImage, 1, 1 );
  
  // Set sliding-specific parameters
  // TODO: Check why there is a massive difference in the gradients between nifty-reg and SuPReMo
  slidingTrafo->SetGapOverlapConstraintWeight( 0.01 / targetImage->nvox );

  slidingTrafo->SetSignedDistanceMapImage( distMapImage );

  // Initialise the current/the only level
  slidingTrafo->InitialiseLevel( 0 );

  // Set the parameters of the sliding transformation need to copy the information from the two individual 
  // transforms first. 
  float* transformParams = (float*)malloc( sizeof( float )*(cpg1Image->nvox + cpg2Image->nvox) );
  memcpy( transformParams, cpg1Image->data, sizeof( float ) * cpg1Image->nvox );
  memcpy( &(transformParams[cpg2Image->nvox]), cpg2Image->data, sizeof( float ) * cpg2Image->nvox );
  slidingTrafo->SetParameters( transformParams, false );

  if (nullptr != expectedDVFImage)
  {
    // First check the DVF, if this is not equal, the transformation will fail
    nifti_image* dvfToTest = slidingTrafo->GetDeformationVectorField( targetImage );

    if (checkImageEquality<SlidingTransformation::PrecisionType>( dvfToTest, expectedDVFImage, EPS_SINGLE, 0, MAX_NUM_DEVIATIONS_ALLOWED ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated DVF image and expected DVF image are not equal." );
      return EXIT_FAILURE;
    }
  }


  // Compare the transformed source image 
  if (nullptr != expectedTransformedSourceImage)
  {
    // Allocate the transformed image
    nifti_image* transformedSourceImage = nifti_copy_nim_info( targetImage );
    transformedSourceImage->data = malloc( sizeof( float ) * targetImage->nvox );
    
    // Use the sliding transformation object to transform the image
    slidingTrafo->TransformImage( sourceImage, transformedSourceImage );

    if (checkImageEquality<SlidingTransformation::PrecisionType>( transformedSourceImage, expectedTransformedSourceImage, EPS_SINGLE, 0, MAX_NUM_DEVIATIONS_ALLOWED ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated transformed source image and expected transformed image are not equal." );
      return EXIT_FAILURE;
    }

    nifti_image_free( transformedSourceImage );
  }

  // Compare the gap/overlap constraint gradient if required
  if (expectedGOCTGrad1Image != nullptr && expectedGOCTGrad2Image != nullptr)
  {
    // Check the number of parameters first
    if (slidingTrafo->GetNumberOfParameters() != (expectedGOCTGrad1Image->nvox + expectedGOCTGrad2Image->nvox))
    {
      supremo_print_error( "Numbers of parameters for gradient comparison not equal between input and generated ones." );
      return EXIT_FAILURE;
    }

    // Generate images for the gradient ()
    nifti_image* constraintGradImage1 = nifti_copy_nim_info( expectedGOCTGrad1Image );
    nifti_image* constraintGradImage2 = nifti_copy_nim_info( expectedGOCTGrad1Image );


    // Get the raw gradient and force it into a newly created image
    SlidingTransformation::PrecisionType* constraintGrad = slidingTrafo->GetConstraintGradientWRTTransformationParameters();
    constraintGradImage1->data = (void*)constraintGrad;
    constraintGradImage2->data = (void*)&(constraintGrad[constraintGradImage1->nvox]);

    reg_io_WriteImageFile( constraintGradImage1, "C:/development/SuPReMo/SuPReMo/motion-test/test-data/xcat/dataGen/supremoGOCTGrad1.nii.gz" );
    reg_io_WriteImageFile( constraintGradImage2, "C:/development/SuPReMo/SuPReMo/motion-test/test-data/xcat/dataGen/supremoGOCTGrad2.nii.gz" );

    if (checkImageEquality<SlidingTransformation::PrecisionType>( constraintGradImage1, expectedGOCTGrad1Image, EPS_SINGLE, 0, 0 ))
    {
      supremo_print_error( "Constraint gradient region 1 unequal." );
      return EXIT_FAILURE;
    }
    if (checkImageEquality<SlidingTransformation::PrecisionType>( constraintGradImage2, expectedGOCTGrad2Image, EPS_SINGLE, 0, 0 ))
    {
      supremo_print_error( "Constraint gradient region 2 unequal." );
      return EXIT_FAILURE;
    }
    
    free( constraintGrad );
  }

  // Free up allocated memory
  free( transformParams );
  
  nifti_image_free( sourceImage );
  nifti_image_free( targetImage );
  nifti_image_free( distMapImage );
  
  if (nullptr != expectedDVFImage)  nifti_image_free( expectedDVFImage );
  if (nullptr != expectedTransformedSourceImage) nifti_image_free( expectedTransformedSourceImage );
  if (nullptr != expectedGOCTGrad1Image) nifti_image_free( expectedGOCTGrad1Image );
  if (nullptr != expectedGOCTGrad2Image) nifti_image_free( expectedGOCTGrad2Image );

  return EXIT_SUCCESS;
}


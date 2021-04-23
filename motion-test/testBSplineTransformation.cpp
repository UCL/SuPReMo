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
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>


#define EPS 0.000001
#define EPS_SINGLE 0.001




int main(int argc, char **argv)
{
  // Print usage information
  if (argc != 7)
  {
    fprintf( stderr, "Usage: %s <cpgImage> <cpgReferenceImage> <referenceImage> <floatingImage> <expectedDVFImage>  <expectedWarpedFloatingImage>\n", argv[0] );
    return EXIT_FAILURE;
  }

  // Get the input parameters
  char *cpgImageName = argv[1];
  char *cpgReferenceImageName = argv[2];
  char *referenceImageName = argv[3];
  char *floatingImageName = argv[4];
  char *expectedDVFImageName = argv[5];
  char *expectedWarpedImageName = argv[6];

  // Read the cpg image 
  nifti_image *cpgImage = reg_io_ReadImageFile( cpgImageName );
  if ( cpgImage == NULL )
  {
    supremo_print_error( "The cpg image could not be read" );
    return EXIT_FAILURE;
  }

  nifti_image *cpgReferenceImage = reg_io_ReadImageFile( cpgReferenceImageName );
  if (cpgReferenceImage == NULL)
  {
    supremo_print_error( "The cpg reference image could not be read" );
    return EXIT_FAILURE;
  }


  // Read the reference image 
  nifti_image *referenceImage = reg_io_ReadImageFile( referenceImageName );
  if ( referenceImage == NULL )
  {
    supremo_print_error( "The reference image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the floating image 
  nifti_image *floatingImage = reg_io_ReadImageFile( floatingImageName );
  if ( floatingImage == NULL )
  {
    supremo_print_error( "The floating image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the expected DVF image 
  nifti_image *expectedDVFImage = reg_io_ReadImageFile( expectedDVFImageName );
  if ( expectedDVFImage == NULL )
  {
    supremo_print_error( "The floating image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the expected warped image 
  nifti_image *expectedWarpedImage = reg_io_ReadImageFile( expectedWarpedImageName );
  if ( expectedWarpedImage == NULL )
  {
    supremo_print_error( "The floating image could not be read" );
    return EXIT_FAILURE;
  }

  // Get the cpg spacing from the image information
  float* cpgSpacing = new float[cpgImage->ndim];

  if ( cpgImage->nu>=2 )
  {
    cpgSpacing[0] = cpgImage->dx;
    cpgSpacing[1] = cpgImage->dy;
  }

  if ( cpgImage->nu > 2 )
  {
    cpgSpacing[2] = cpgImage->dz;
  }

  // Generate the b-spline transformation object
  std::shared_ptr<BSplineTransformation> bsplTrafo = std::make_shared<BSplineTransformation>( cpgReferenceImage, 1, cpgSpacing );
   
  // Check that the generated B-spline transformation and the input cpg have the same number of parameters
  if ( bsplTrafo->GetNumberOfParameters() != cpgImage->nvox )
  {
    supremo_print_error( "Wrong number of parameters in BSpline-Transformation" );
    return EXIT_FAILURE;
  }

  // Since the data is copied over using memcpy 
  // set the parametres from the input controlpoint grid
  bsplTrafo->SetParameters( (float*) cpgImage->data, false );

  // First check the DVF, if this is not equal, the transformation will fail
  nifti_image* dvfToTest = bsplTrafo->GetDeformationVectorField( referenceImage );

  if ( checkImageEquality<BSplineTransformation::PrecisionType>( dvfToTest, expectedDVFImage, EPS_SINGLE ) == EXIT_FAILURE )
  {
    supremo_print_error( "Generated DVF image and expected DVF image are not equal." );
    return EXIT_FAILURE;
  }

  // Check the transformed image for equality with the expected outcome
  nifti_image* transformedImg;
  transformedImg = nifti_copy_nim_info( referenceImage );
  transformedImg->data = malloc( transformedImg->nvox * transformedImg->nbyper );
  bsplTrafo->TransformImage( floatingImage, transformedImg );

  if ( checkImageEquality<BSplineTransformation::PrecisionType>( transformedImg, expectedWarpedImage, EPS_SINGLE ) == EXIT_FAILURE )
  {
    supremo_print_error( "Generated DVF image and expected DVF image are not equal." );
    return EXIT_FAILURE;
  }

  
  //// Check the inversion of the transformation
  //// Go from BSpline grid to DVF and back again
  //// Use the DVF as above and convert it back to CPGs
  //  
  //// Generate an image to check the equality.
  //reg_getDisplacementFromDeformation( dvfToTest );
  //reg_getDisplacementFromDeformation( cpgImage );
  //nifti_image* cpgImgRecovered = nifti_copy_nim_info( cpgImage );
  //cpgImgRecovered->data = (void*)bsplTrafo->GetDVFGradientWRTTransformationParameters( dvfToTest, referenceImage );
  //
  //reg_io_WriteImageFile( cpgImgRecovered, "mycpgImgRecovered.nii.gz" );
  //reg_io_WriteImageFile( dvfToTest, "mydvfToTest.nii.gz" );
  //reg_io_WriteImageFile( cpgImage, "mycpgImage.nii.gz" );
  //if (checkImageEquality<BSplineTransformation::PrecisionType>( cpgImgRecovered, cpgImage, EPS_SINGLE )==EXIT_FAILURE)
  //{
  //  nmm_print_error( "Generated DVF image and expected DVF image are not equal." );
  //  return EXIT_FAILURE;
  //}
  //
  
  
  
  
  
  //// Test the adjoint transformation
  //// Show that the following property holds:
  //// < Tx, y > = < x, T* y >
  //// Tx == transformedImg
  //// x  == floatingImage
  //// y  == 




  //// Allocate the weight images
  //nifti_image* sourceWeightsImg = nifti_copy_nim_info( transformedImg );
  //sourceWeightsImg->data = calloc( sourceWeightsImg->nvox, sourceWeightsImg->nbyper );

  //nifti_image* adjointWeightsImg = nifti_copy_nim_info( transformedImg );
  //adjointWeightsImg->data = calloc( sourceWeightsImg->nvox, sourceWeightsImg->nbyper );

  //// Fill the source weights with 1
  //float* sourceWeightsPointer = static_cast<float*>(sourceWeightsImg->data);
  //for (unsigned int uiI = 0; uiI < sourceWeightsImg->nvox; ++uiI)
  //{
  //  sourceWeightsPointer[uiI] = 1.f;
  //}

  //// allocate a result image 
  //nifti_image* adjointTransformedImg = nifti_copy_nim_info( transformedImg );
  //adjointTransformedImg->data = calloc( adjointTransformedImg->nvox, adjointTransformedImg->nbyper );

  //// Perform the adjoint transformation 
  //bsplTrafo->TransformImageAdjoint( transformedImg, sourceWeightsImg, adjointTransformedImg, adjointWeightsImg );
  //
  //// Todo: Normalisation with adjoint weights
  //float* adjWeigthsPtr = static_cast<float*>(adjointWeightsImg->data);
  //float* adjImgPtr = static_cast<float*>(adjointTransformedImg->data);
  //
  //for (unsigned int uiI = 0; uiI < sourceWeightsImg->nvox; ++uiI)
  //{
  //  if (adjImgPtr[uiI] != 0)
  //    adjImgPtr[uiI] = adjImgPtr[uiI] / adjWeigthsPtr[uiI];
  //  else
  //    adjImgPtr[uiI] = 0.f;
  //}

  ///// ToDo: Implement test according to https://math.stackexchange.com/questions/262930/adjoint-operators-and-inner-product-spaces?answertab=votes#tab-top
  ///// As per discussion with Richardk on 8th Jan 2020.

  //reg_io_WriteImageFile( adjointTransformedImg, "C:/data/testAdjointImgOut.nii.gz" );
  //reg_io_WriteImageFile( transformedImg, "C:/data/testTransformedImgOut.nii.gz" );
  return EXIT_SUCCESS;
}

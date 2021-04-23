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





//--------------------------------------
// Includes
//--------------------------------------
#include "SSDImageSimilarity.h"
#include "Supremo.h"




//----------------------------------------
// SSDImageSimilarity::SSDImageSimilarity
//----------------------------------------
SSDImageSimilarity::SSDImageSimilarity()
{}




//-----------------------------------------
// SSDImageSimilarity::~SSDImageSimilarity
//-----------------------------------------
SSDImageSimilarity::~SSDImageSimilarity()
{}




//--------------------------------------------------------
// SSDImageSimilarity::GetSimilarityMeasureValueForImages
//--------------------------------------------------------
double SSDImageSimilarity::GetSimilarityMeasureValueForImages( nifti_image* img1In,
                                                               nifti_image* img2In )
{
  /// \todo Perform checks that images are of equal size, and equal data type (equal geometry?)

  // Generate the values requried for using nifti_reg function
  // Only consider the first time point in the image
  //double* activeTimePoints = new double[img1In->nt]();
  //activeTimePoints[0] = 1.0;

  // fake mask all set to zero means all voxels are being used
  //int* fakeMask = new int[img1In->nvox]();

  // value not documented in niftyReg.
  //float* currentValue = new float[img1In->nt]();
  double ssdValToReturn = GetSSDValue( img1In, img2In );
  //double ssdValToReturn = reg_getSSDValue<PrecisionType>( img1In, img2In, activeTimePoints, nullptr, fakeMask, currentValue, nullptr );
  
  // Clean up
  //delete[] activeTimePoints;
  //delete[] fakeMask;
  //delete[] currentValue;

  // And exit  
  return ssdValToReturn;
}




//----------------------------------------------------
// SSDImageSimilarity::GetSimilarityGradientWRTVoxels
//----------------------------------------------------
void SSDImageSimilarity::GetSimilarityGradientWRTVoxels( nifti_image* referenceImg,
														 nifti_image* sourceImg,
														 nifti_image* similarityGradientWRTVoxelsOutputImage)
{

	// Check image if images are of same data type
	if (sourceImg->datatype != referenceImg->datatype)
	{
		supremo_print_error("Reference and souce image have to have the same data type");
		supremo_exit(EXIT_FAILURE, __FILE__, __LINE__);
	}

	// Generate the output image
	/// \todo Decide how to deal with multiple time points in images
	PrecisionType* refDataPtr = (PrecisionType*)referenceImg->data;
	PrecisionType* srcDataPtr = (PrecisionType*)sourceImg->data;
	PrecisionType* gradDataPtr = (PrecisionType*)similarityGradientWRTVoxelsOutputImage->data;

	for (long int voxel = 0; voxel < referenceImg->nvox; ++voxel)
	{
		PrecisionType referenceValue, sourceValuel;
		referenceValue = refDataPtr[voxel] * referenceImg->scl_slope + referenceImg->scl_inter;
		sourceValuel = srcDataPtr[voxel] * sourceImg->scl_slope + sourceImg->scl_inter;
		gradDataPtr[voxel] = -2.0f * (referenceValue - sourceValuel);
	}

	return;
}





double SSDImageSimilarity::GetSSDValue( nifti_image * referenceImage, nifti_image * sourceImage )
{
  int voxel;
  int voxelNumber = (int)referenceImage->nx * referenceImage->ny * referenceImage->nz;

  // Create pointers to the reference and warped image data
  PrecisionType *referencePtr = static_cast<PrecisionType *>(referenceImage->data);
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);

  double SSD_global = 0.0, n = 0.0;
  double targetValue, resultValue;

  // Loop over the different time points
  double SSD_local = 0.;

#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
    shared(referenceImage, referencePtr, sourcePtr, \
          voxelNumber) \
    private(voxel, targetValue, resultValue) \
reduction(+:SSD_local) \
reduction(+:n)
#endif
  for (voxel = 0; voxel < voxelNumber; ++voxel)
  {
    // Ensure that both reference and source values are defined
    targetValue = (double)(referencePtr[voxel] * referenceImage->scl_slope + referenceImage->scl_inter);
    resultValue = (double)(sourcePtr[voxel] * referenceImage->scl_slope + referenceImage->scl_inter);
    if (targetValue == targetValue && resultValue == resultValue)
    {
      double diff = targetValue - resultValue;
      SSD_local += diff * diff;
      n += 1.0;
    }
  }
  SSD_global -= SSD_local / n;
  
  return SSD_global;
}

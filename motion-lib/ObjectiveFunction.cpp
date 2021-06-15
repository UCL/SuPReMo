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





//----------
// Includes
//----------
#include "ObjectiveFunction.h"
#include "_reg_ReadWriteImage.h"


//--------------------------------------
// ObjectiveFunction::ObjectiveFunction
//--------------------------------------
ObjectiveFunction::ObjectiveFunction() :
  correspondenceModel( nullptr ),
  imageSimilarity( nullptr ),
  referenceStateImage( nullptr ),
  dynamicDataType( SAME_RES_AS_STATIC ),
  imageAcquisition( nullptr )
{}




//---------------------------------------
// ObjectiveFunction::~ObjectiveFunction
//---------------------------------------
ObjectiveFunction::~ObjectiveFunction()
{}




//-----------------------------------
// ObjectiveFunction::GetMaxStepSize
//-----------------------------------
ObjectiveFunction::PrecisionType ObjectiveFunction::GetMaxStepSize()
{
  // Check that the surrogate signal was set
  if (this->surrogateSignals.empty())
  {
    supremo_print_error( "Surrogate signal must not be empty when calling ObjectiveFunction::GetMaxStepSize()." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Check that the reference state image was set
  if (nullptr == this->referenceStateImage)
  {
    supremo_print_error( "Reference state image must be set before calling ObjectiveFunction::GetMaxStepSize()." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Calculate the maximum voxel size 
  PrecisionType maxStepSize = this->referenceStateImage->dx > this->referenceStateImage->dy ? this->referenceStateImage->dx : this->referenceStateImage->dy;
  if (this->referenceStateImage->ndim > 2)
  {
    maxStepSize = maxStepSize > this->referenceStateImage->dz ? maxStepSize : this->referenceStateImage->dz;
  }

  // Calculate the maximum absolute surrogate signal
  PrecisionType maxSurrSignalValue = 0.;

  for (size_t i = 0; i < this->surrogateSignals.size(); ++i)
  {
    for (size_t j = 0; j < this->surrogateSignals[i].size(); j++) 
    {
      maxSurrSignalValue = maxSurrSignalValue > fabs( this->surrogateSignals[i][j] ) ? maxSurrSignalValue : fabs( this->surrogateSignals[i][j] );
    }
  }

  maxStepSize = maxStepSize / maxSurrSignalValue;

  return maxStepSize;
}




//-------------------------------------------
// ObjectiveFunction::SetCorrespondenceModel
//-------------------------------------------
void ObjectiveFunction::SetCorrespondenceModel( const std::shared_ptr<CorrespondenceModel> & correspondenceModelIn )
{
  // Copy the shared pointer of the correspondence model
  this->correspondenceModel = correspondenceModelIn;
  
  // update the similarity weight based on the transformation contained within the correspondence model
  this->similarityWeight = (PrecisionType) 1. - this->correspondenceModel->GetTransformation()->GetSumOfPenaltyWeights();
}




//-------------------------------------------
// ObjectiveFunction::SetReferenceStateImage
//-------------------------------------------
void ObjectiveFunction::SetReferenceStateImage( nifti_image* referenceSateImageIn )
{
  this->referenceStateImage = referenceSateImageIn;
}




//-------------------------------------------
// ObjectiveFunction::SetDynamicImages
//-------------------------------------------
void ObjectiveFunction::SetDynamicImages( const std::vector<nifti_image*>& dynamicImagesIn, t_dynamicData dynamicDataTypeIn )
{
  this->dynamicImages = dynamicImagesIn;
  this->dynamicDataType = dynamicDataTypeIn;
}




//----------------------------------------
// ObjectiveFunction::SetImageAcquisition
//----------------------------------------
void ObjectiveFunction::SetImageAcquisition( const std::shared_ptr<ImageAcquisition>& imageAcquisitionIn )
{
  // Simply set the shared image acquisition pointer to the member variable.
  this->imageAcquisition = imageAcquisitionIn;
}




//----------------------------------------
// ObjectiveFunction::SetSurrogateSignals
//----------------------------------------
void ObjectiveFunction::SetSurrogateSignals( const ObjectiveFunction::SurrogateSignalType & surrogateSignalsIn )
{
  this->surrogateSignals = surrogateSignalsIn;
}




//----------------------------------------------
// ObjectiveFunction::SetImageSimilarityMeasure
//----------------------------------------------
void ObjectiveFunction::SetSimilarityMeasure( const std::shared_ptr<ImageSimilarity>& imageSimilarityIn )
{
  this->imageSimilarity = imageSimilarityIn;
}




//-----------------------------
// ObjectiveFunction::GetValue
//-----------------------------
ObjectiveFunction::PrecisionType ObjectiveFunction::GetValue( const ObjectiveFunction::PrecisionType * parametersIn )
{
  // Check that the number of time points of the surrogate signals corresponds those of the dynamic images
  if ( this->surrogateSignals.size() != this->dynamicImages.size() )
  {
    supremo_print_error("Number of time points of dyanmic images and surrogate signals do not match.");
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Sum of the constraint values
  double currentConstraintSum = 0.0;
  
  // Sum of the weighted similarity measures
  double currentWeightedSimilarityMeasure = 0.0;
  
  // Set the correspondence model to the correct parameters
  this->correspondenceModel->SetParameters( parametersIn );

  // Iterate over all time points
  for (int iTimePoint = 0; iTimePoint < this->surrogateSignals.size(); ++iTimePoint)
  {
    // Get the current transformation (based on the surrogate signal) and calculate the 
    std::shared_ptr<Transformation> curTransformation = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[iTimePoint] );

    // The image acquisition needs to determine, what is the minimum size of the warped reference
    // state image to warp only as much of the reference state image actually necessary.
    // Note other image acquisition methods might need to provide additional information
    nifti_image* warpedRefStateImage = this->imageAcquisition->AllocateMinimumSizeImgInFullImgSpace( this->referenceStateImage, this->dynamicImages[iTimePoint] );
    
    // transform the image, the result will be written to warpedRefStateImage
    curTransformation->TransformImage( this->referenceStateImage, warpedRefStateImage );

    // simulate the image acquisition
    nifti_image* simulatedDynamicImage = this->imageAcquisition->SimulateImageAcquisition( warpedRefStateImage, this->dynamicImages[iTimePoint] );

    // Calculate the similarity value between the simulated dynamic image and the dynamic image
    // Note that the weighting takes place here. Weighting was set when correspondence model and 
    // contained transformation were set.
    double similarityValue = imageSimilarity->GetSimilarityMeasureValueForImages( simulatedDynamicImage, this->dynamicImages[iTimePoint] );
    if ( similarityValue == similarityValue )
    {
      currentWeightedSimilarityMeasure += this->similarityWeight * similarityValue;
    } 

    // Calculate the constraint value
    currentConstraintSum += curTransformation->GetConstraintValue();

    // Clean up the reference state image and the simulated dynamic image
    if (nullptr != warpedRefStateImage)   nifti_image_free( warpedRefStateImage );
    if ((nullptr != simulatedDynamicImage) && (warpedRefStateImage != simulatedDynamicImage)) nifti_image_free( simulatedDynamicImage );
  }

  return currentWeightedSimilarityMeasure - currentConstraintSum;
}




//--------------------------------
// ObjectiveFunction::GetGradient
//--------------------------------
void ObjectiveFunction::GetGradient( const PrecisionType* parametersIn, PrecisionType* gradientOut, bool normaliseGradient )
{
  // Set the correspondence model to the correct parameters
  this->correspondenceModel->SetParameters( parametersIn );

  // Reset the input gradient to zero
  for (unsigned int iDOF = 0; iDOF < this->GetNumberOfParameters(); ++iDOF)
  {
    gradientOut[iDOF] = (PrecisionType) 0.0;
  }
  
  // The gradient of the respiratory motion model is calculated by summing over all time points
  // Iterate over all time points
  for ( int iTimePoint = 0; iTimePoint < this->surrogateSignals.size(); ++iTimePoint )
  {
    // Get the current transformation (based on the surrogate signal) and calculate the 
    // constraint value
    std::shared_ptr<Transformation> curTransformation = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[iTimePoint] );

    // Prepare the two images involved at this time point
    nifti_image* curDynamicImage = this->dynamicImages[iTimePoint];

    // The image acquisition needs to determine, what is the minimum size of the warped reference
    // state image to warp only as much of the reference state image actually necessary.
    // Note other image acquisition methods might need to provide additional information
    nifti_image* warpedRefStateImage = this->imageAcquisition->AllocateMinimumSizeImgInFullImgSpace( this->referenceStateImage, this->dynamicImages[iTimePoint] );

    // transform the image, the result will be written to warpedRefStateImage
    curTransformation->TransformImage( this->referenceStateImage, warpedRefStateImage );

    // simulate the image acquisition
    nifti_image* simulatedDynamicImage = this->imageAcquisition->SimulateImageAcquisition( warpedRefStateImage, this->dynamicImages[iTimePoint] );

    // Allocate output image...
    nifti_image* simGradImgWRTSimDynVoxels = nifti_copy_nim_info( curDynamicImage );
    simGradImgWRTSimDynVoxels->data = (void*)calloc( simGradImgWRTSimDynVoxels->nvox, simGradImgWRTSimDynVoxels->nbyper );
    // ... and calculate the image similarity gradient wrt. the voxel intensities
    this->imageSimilarity->GetSimilarityGradientWRTVoxels( curDynamicImage, simulatedDynamicImage, simGradImgWRTSimDynVoxels );

    // use adjoint of image acquisition to transform measure gradient into warped image space
    // store in dynamicAfterAdjointImage
    this->imageAcquisition->CalculateAdjoint( this->referenceStateImage, simGradImgWRTSimDynVoxels );
    nifti_image* dynamicAfterAdjointImage = this->imageAcquisition->GetImageAfterAdjoint();

    // Allocate output image
    nifti_image* imageGradientWRTDVF = nifti_copy_nim_info( curTransformation->GetDeformationVectorField( dynamicAfterAdjointImage ) );
    imageGradientWRTDVF->data = calloc( imageGradientWRTDVF->nvox, imageGradientWRTDVF->nbyper );

    curTransformation->GetImageGradientWRTDVF( this->referenceStateImage, imageGradientWRTDVF );

    // Allocate "voxelBasedMeasureGrad"
    nifti_image* imageGradientWRTDVFTimesAdjoint = nifti_copy_nim_info( imageGradientWRTDVF );
    imageGradientWRTDVFTimesAdjoint->data = (void *) calloc( imageGradientWRTDVFTimesAdjoint->nvox,
                                                             imageGradientWRTDVFTimesAdjoint->nbyper );

    //then it is multiplied by the measure gradient
    long voxel;
    long numVox = dynamicAfterAdjointImage->nx * dynamicAfterAdjointImage->ny * dynamicAfterAdjointImage->nz;
    PrecisionType *dynamicAfterAdjointDataPtr = static_cast<PrecisionType *>(dynamicAfterAdjointImage->data);
    PrecisionType *imgGradWRTDVFDataPtrX = static_cast<PrecisionType *>(imageGradientWRTDVF->data);
    PrecisionType *imgGradWRTDVFDataPtrY = &imgGradWRTDVFDataPtrX[numVox];
    PrecisionType *imgGradWRTDVFDataPtrZ = nullptr;
    if ( imageGradientWRTDVF->nu == 3 )
      imgGradWRTDVFDataPtrZ = &imgGradWRTDVFDataPtrY[numVox];
    PrecisionType *imgGradWRTDVFTimesAdjDataPtrX = static_cast<PrecisionType *>(imageGradientWRTDVFTimesAdjoint->data);
    PrecisionType *imgGradWRTDVFTimesAdjDataPtrY = &imgGradWRTDVFTimesAdjDataPtrX[numVox];
    PrecisionType *imgGradWRTDVFTimesAdjDataPtrZ = nullptr;
    if ( imageGradientWRTDVFTimesAdjoint->nu == 3 )
      imgGradWRTDVFTimesAdjDataPtrZ = &imgGradWRTDVFTimesAdjDataPtrY[numVox];
    PrecisionType tmpVal;
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
         shared(dynamicAfterAdjointDataPtr, imgGradWRTDVFDataPtrX, imgGradWRTDVFDataPtrY, imgGradWRTDVFDataPtrZ, \
			imgGradWRTDVFTimesAdjDataPtrX, imgGradWRTDVFTimesAdjDataPtrY, imgGradWRTDVFTimesAdjDataPtrZ, numVox) \
		 private(voxel, tmpVal)
#endif
    for ( voxel = 0; voxel < numVox; voxel++ )
    {
      tmpVal = dynamicAfterAdjointDataPtr[voxel] * imgGradWRTDVFDataPtrX[voxel];
      if ( !isnan(tmpVal) ) 
        imgGradWRTDVFTimesAdjDataPtrX[voxel] = tmpVal;
      tmpVal = dynamicAfterAdjointDataPtr[voxel] * imgGradWRTDVFDataPtrY[voxel];
      if ( !isnan(tmpVal) )
        imgGradWRTDVFTimesAdjDataPtrY[voxel] = tmpVal;
      if (imgGradWRTDVFDataPtrZ != nullptr)
      {
        tmpVal = dynamicAfterAdjointDataPtr[voxel] * imgGradWRTDVFDataPtrZ[voxel];
        if (!isnan( tmpVal ))
          imgGradWRTDVFTimesAdjDataPtrZ[voxel] = tmpVal;
      }
    }

    // Calculate the transformation gradient 
    PrecisionType* transformationGradient = curTransformation->GetDVFGradientWRTTransformationParameters( imageGradientWRTDVFTimesAdjoint );

    // Get the gradient in the direction of the selected constraints
    // Constraint weights are applied by the transformation
    PrecisionType* constraintGradient = curTransformation->GetConstraintGradientWRTTransformationParameters();

    // Multiply the transformation gradient by the similarity weight and add the constraint gradient
    unsigned int totalNumberOfTrafoParameters = curTransformation->GetNumberOfParameters();
    for (unsigned int nParam = 0; nParam < totalNumberOfTrafoParameters; ++nParam)
    {
      transformationGradient[nParam] = transformationGradient[nParam] * this->similarityWeight + constraintGradient[nParam];
    }

    // Calculate the gradient update and add it to the gradientOut
    this->correspondenceModel->GetTransformationParameterGradientWRTModelParameters( transformationGradient, this->surrogateSignals[iTimePoint], gradientOut );

    // Clean up for this time-point
    nifti_image_free( simGradImgWRTSimDynVoxels );
    
    // If no image acquisition was simulated, these will be the same, so need to free up only one and set the other one to null. 
    if (warpedRefStateImage == simulatedDynamicImage)
    {
      nifti_image_free( warpedRefStateImage );
      warpedRefStateImage = nullptr;
      simulatedDynamicImage = nullptr;
    }
    else
    {
      nifti_image_free( warpedRefStateImage );
      nifti_image_free( simulatedDynamicImage );
      warpedRefStateImage = nullptr;
      simulatedDynamicImage = nullptr;
    }
    nifti_image_free( imageGradientWRTDVF );
    nifti_image_free( imageGradientWRTDVFTimesAdjoint );

    free( transformationGradient );
    free( constraintGradient );
  }  

  // Normalise the gradient if required
  if (normaliseGradient)
  {
    PrecisionType gradLength = this->correspondenceModel->GetMaxParameterLength( gradientOut );
#ifdef _DEBUG
    std::cout << "Normalising gradient to maximum length of " << gradLength << std::endl;
#endif
    for (unsigned int iDOF = 0; iDOF < this->GetNumberOfParameters(); ++iDOF)
    {
      gradientOut[iDOF] = gradientOut[iDOF] / gradLength;
    }
  }
}



//---------------------------------------
// ObjectiveFunction::GetGradientAsImage
//---------------------------------------
std::vector<nifti_image*> ObjectiveFunction::GetGradientAsImage( const PrecisionType* parametersIn, bool normaliseGradient )
{
  // Overview: 
  // The objective function gradient has the same size as the correspondence model. Hence, the 
  // correspondence model functionality is used to generate an image first, then replace the contents
  // with the calculated gradient.
  
  // Generate the output image from the correspondence model
  // This will have the same dimension and data type as the gradient
  std::vector<nifti_image*> outGradientImages = this->correspondenceModel->GetCorrespondenceModelAsImage();
  
  // Since an unknown number of images can be returned, we have to split the data up in the order the images are returned
   
  // The total number of parameters must coincide with the number of voxels in the images returned by the correspondence model
  unsigned int nParameters = this->GetNumberOfParameters();
  unsigned int totalNumOfVoxels = 0;

  for (size_t nImg = 0; nImg < outGradientImages.size(); nImg++)
  {
    totalNumOfVoxels += outGradientImages[nImg]->nvox;
  }

  if (nParameters != totalNumOfVoxels)
  {
    supremo_print_error( "Number of voxels in output image is not equal to number of parameters of the gradient." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Save the gradient to an intermediate array
  PrecisionType* gradData = static_cast<PrecisionType*>(malloc( sizeof( PrecisionType ) * nParameters ));
  this->GetGradient( parametersIn, gradData, normaliseGradient );

  // From there copy it over to the images provided by the correspondence model
  unsigned int curVoxStartIdx = 0;
  
  for (size_t nImg = 0; nImg < outGradientImages.size(); nImg++)
  {
    memcpy( outGradientImages[nImg]->data, &(gradData[curVoxStartIdx]), outGradientImages[nImg]->nvox * sizeof( PrecisionType ) );
      curVoxStartIdx += outGradientImages[nImg]->nvox;
  }

  // Clear the memory
  free( gradData );
  
  return outGradientImages;
}




//------------------------------------------
// ObjectiveFunction::GetNumberOfParameters
//------------------------------------------
unsigned int ObjectiveFunction::GetNumberOfParameters()
{
	return this->correspondenceModel->GetNumberOfParameters();
}


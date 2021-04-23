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
// includes
//----------
#include "MoCoReconWeightedAverage.h"
#include "Supremo.h"



//----------------------------------------------------
// MoCoReconWeightedAverage::MoCoReconWeightedAverage
//----------------------------------------------------
MoCoReconWeightedAverage::MoCoReconWeightedAverage() :
  weightsImageForReconstruction( nullptr )
{
  // Indicate the behaviour of the average weighted reconstruction.  
  // Repeatedly calling with same correspondence model will not change the result of the MCR.
  this->repeatedUpdateChangesMCRImage = false;
}




//-----------------------------------------------------
// MoCoReconWeightedAverage::~MoCoReconWeightedAverage
//-----------------------------------------------------
MoCoReconWeightedAverage::~MoCoReconWeightedAverage()
{
  this->ClearWeightsImageForReconstruction();
  return;
}




//----------------------------------------------------
// MoCoReconWeightedAverage::Update
//----------------------------------------------------
void MoCoReconWeightedAverage::Update()
{
  std::cout << "Performing motion-compensated reconstruction (average weighting)" << std::endl;

  // Need to allocate a weights image and reconstructed image that both have the size 
  // of the image defining the reconstruction geometry.
  // Do not need to set images explicitly to zero, since allocation takes care of this
  this->AllocateReconstructedImage();
  this->AllocateWeightsImageForReconstruction();
  
  for (size_t i = 0; i < this->numberOfDynamicImages; ++i) 
  {
    // Get the transformation from the surrogate signal
    std::shared_ptr<Transformation> curTrafo = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[i] );

    // Get the adjoint of the acquired image. 
    this->imageAcquisition->CalculateAdjoint( this->reconstructionGeometryImage, this->dynamicImages[i] );

    // Use the implementation of Transformation class to calculate the adjoint/push transformation
    curTrafo->TransformImageAdjoint( 
      this->imageAcquisition->GetImageAfterAdjoint(), 
      this->imageAcquisition->GetWeightsImageAfterAdjoint(),
      this->reconstructedImage, 
      this->weightsImageForReconstruction );
  }

  // divide reconstructed image by weights image to get weighted average
  // any voxels with 0 weight should be set to padding value
  PrecisionType *recoPtr = static_cast<PrecisionType *>(this->reconstructedImage->data);
  PrecisionType *weightPtr = static_cast<PrecisionType *>(this->weightsImageForReconstruction->data);
  
  // Get the padding value from the transformaiton used
  Transformation::PrecisionType paddingValue = this->correspondenceModel->GetTransformation()->GetPaddingValue();
  
  for (int v = 0; v < this->reconstructedImage->nvox; ++v)
  {
    // Only normalise by the weights image if unequal to zero
    if (weightPtr[v] == 0)
      recoPtr[v] = paddingValue;
    else
      recoPtr[v] = recoPtr[v] / weightPtr[v];
  }

  return;
}




//-----------------------------------------------------------------
// MoCoReconWeightedAverage::AllocateWeightsImageForReconstruction
//-----------------------------------------------------------------
void MoCoReconWeightedAverage::AllocateWeightsImageForReconstruction()
{
  // Check that the image that defines the geometry exists
  if (nullptr == this->reconstructionGeometryImage) 
  {
    supremo_print_error( "Geometry image has to be set before allocating weights image." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Clear the weights image
  this->ClearWeightsImageForReconstruction();

  // Generate a new one with precision type, and scl-slope and intersection of 1.0 and 0.0 respectively
  // Note: double precision is not tested 
  this->weightsImageForReconstruction = nifti_copy_nim_info( this->reconstructionGeometryImage );
  this->weightsImageForReconstruction->scl_slope = 1.0f;
  this->weightsImageForReconstruction->scl_inter = 0.0f;
  this->weightsImageForReconstruction->nbyper = sizeof( PrecisionType );
  
  if (sizeof( PrecisionType ) == sizeof( float ))
    this->weightsImageForReconstruction->datatype = NIFTI_TYPE_FLOAT32;
  else
    this->weightsImageForReconstruction->datatype = NIFTI_TYPE_FLOAT64;

  this->weightsImageForReconstruction->data = (void *)calloc( this->weightsImageForReconstruction->nvox, this->weightsImageForReconstruction->nbyper );
}




//--------------------------------------------------------------
// MoCoReconWeightedAverage::ClearWeightsImageForReconstruction
//--------------------------------------------------------------
void MoCoReconWeightedAverage::ClearWeightsImageForReconstruction()
{
  // Only clear the weights image if it was not empty before
  if (nullptr != this->weightsImageForReconstruction)
  {
    nifti_image_free( this->weightsImageForReconstruction );
    this->weightsImageForReconstruction = nullptr;
  }
}

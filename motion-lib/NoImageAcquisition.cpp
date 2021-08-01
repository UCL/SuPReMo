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





#include "NoImageAcquisition.h"


//----------------------------------------
// NoImageAcquisition::NoImageAcquisition
//----------------------------------------
NoImageAcquisition::NoImageAcquisition()
{}




//-----------------------------------------
// NoImageAcquisition::~NoImageAcquisition
//-----------------------------------------
NoImageAcquisition::~NoImageAcquisition()
{}




//----------------------------------------------
// NoImageAcquisition::SimulateImageAcquisition
//----------------------------------------------
nifti_image* NoImageAcquisition::SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace, unsigned int dynamicImageTimePoint )
{
  // Since no image acqusition is simulated, the image in full space is simply returned. 
  return imgInFullImgSpace;
}




//--------------------------------------
// NoImageAcquisition::CalculateAdjoint
//--------------------------------------
void NoImageAcquisition::CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint )
{
  // The current image in acquisition and full image space as well as the dynamic image time point is remembered for the allocation of the image after adjoint
  this->curImageInAcquisitionSpace = imgInAcquisitionSpace;
  this->curImageInFullImgSpace = imgInFullImgSpace;
  this->curDynamicImageTimePoint = dynamicImageTimePoint;
  
  // Allocate the image after applying the adjoint and the corresponding weights image
  this->AllocateImageAfterAdjoint();
  this->AllocateWeightsImageAfterAdjoint();

  // when not simulating any image acquisition, just copy current image in acquisition space into the image after applying the adjoint
  memcpy( this->imageAfterAdjoint->data, this->curImageInAcquisitionSpace->data, this->imageAfterAdjoint->nvox * this->imageAfterAdjoint->nbyper );
  
  // all voxels in the weights image after adjoint are set to 1
  PrecisionType* tmpPtr = static_cast<PrecisionType *>(this->weightsImageAfterAdjoint->data);
  for (int v = 0; v < this->weightsImageAfterAdjoint->nvox; v++)
  {
    tmpPtr[v] = 1.0;
  }
}




//----------------------------------------------------------
// NoImageAcquisition::AllocateMinimumSizeImgInFullImgSpace
//----------------------------------------------------------
nifti_image * NoImageAcquisition::AllocateMinimumSizeImgInFullImgSpace( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace, unsigned int dynamicImageTimePoint )
{
  // If no image acquisition is simulated, the allocated image has to have the size of the image in acquisition space.  
  nifti_image* minSizedImgInFullImgSpace;

  // Fill in the header information 
  minSizedImgInFullImgSpace = nifti_copy_nim_info( imgInAcquisitionSpace );
  
  // and allocate the image data memory
  minSizedImgInFullImgSpace->data = (void *)calloc( minSizedImgInFullImgSpace->nvox, minSizedImgInFullImgSpace->nbyper );
  return minSizedImgInFullImgSpace;
}




//-----------------------------------------------
// NoImageAcquisition::AllocateImageAfterAdjoint
//-----------------------------------------------
void NoImageAcquisition::AllocateImageAfterAdjoint()
{
  // Delete the previous image after applying the adjoint of the image acquisition
  this->ClearImageAfterAdjoint();
  
  // Allocate an image of the same size and geometry as the acquired image
  // Use the implemented, exposed functionality to allocate the image after adjoint to avoid code duplication
  this->imageAfterAdjoint = this->AllocateMinimumSizeImgInFullImgSpace( this->curImageInFullImgSpace, this->curImageInAcquisitionSpace, this->curDynamicImageTimePoint );

  /// Note: This implementation is different to the reference implementaiton in reg_resp which adjusted the image size of the adjoint 
  ///       regardless of which acquisition simulation was selected. This could lead to a segmentation fault in some specific cases.
  ///       As a result, here the image size is not adjusted for the no-image-acquisition class. 
}


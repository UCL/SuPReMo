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




#include "MoCoReconSuperResolutionIBP.h"
#include <iostream>
#include "Supremo.h"


MoCoReconSuperResolutionIBP::MoCoReconSuperResolutionIBP( unsigned int maxIterationNumberIn, bool updateReconstruction ) :
  maxIterationNumber( maxIterationNumberIn ),
  updateReconstruction( updateReconstruction ),
  correctionImage( nullptr ),
  weightsImageForReconstruction( nullptr )
{
  // If we update the previous result, (updateReconstruction==true) then a repated call of Update() will change the resulting MCR.
  this->repeatedUpdateChangesMCRImage = updateReconstruction;

  return;
}





MoCoReconSuperResolutionIBP::~MoCoReconSuperResolutionIBP()
{
  this->ClearCorrectionImage();
  this->ClearReconstructedImage();
  this->ClearWeightsImageForReconstruction();
}





//-------------------------------------
// MoCoReconSuperResolutionIBP::Update
//-------------------------------------
void MoCoReconSuperResolutionIBP::Update()
{
  // Print out basic information about the motion-compensated reconstruction
  std::cout << "Performing motion-compensated image reconstruction." << std::endl
    << "using iterative back-projection super-resolution method..." << std::endl;
  
  /// Todo:  Need to check that the image acquisition is correct. Maybe best to do this extenrally, since this interdependency is not helpful
  
  // Only allocate a new reconstructed image if we do not want to update the previous reconstruction
  // When restarting the reconstruction, need to always allocate a new reconstructed image
  // when updating, then we only need to allocate a new image, if the old one cannot be updated. 
  if (this->updateReconstruction)
  {
    if (!this->CheckReconstrutedImageCanBeUpdated()) this->AllocateReconstructedImage();
  }
  else
  {
    // restarting reconstruction, always allocate a new reconstructed image, discard previous result
    this->AllocateReconstructedImage();
  }

  // Allocate the correction image and the 
  this->AllocateCorrectionImage();
  this->AllocateWeightsImageForReconstruction();

  // Tracking the SSD value over the iterations. Initialised value of -1 will be used 
  // further down to indicate the first evaluation. 
  PrecisionType prevSSD = -1.;

  // Need an image that stores the resulting reconstruction from the previous iteration. 
  // Initialise with same header info as reconstruction geometry iamge (recon) image and allocate space for data
  nifti_image* prevReconImage = nifti_copy_nim_info( this->reconstructionGeometryImage );
  prevReconImage->data = calloc( prevReconImage->nvox, prevReconImage->nbyper );

  PrecisionType paddingVal;

  // perform up to max iterations - note, start at -1 as correction from last iteration is
  // always removed afterwards, so first iteration effectively leaves recon image unchanged
  int iter;
  for (iter = -1; iter < static_cast<int>(this->maxIterationNumber); ++iter)
  {
    // Set values in weight image and correction image to 0
    PrecisionType *corPtr = static_cast<PrecisionType *>(this->correctionImage->data);
    PrecisionType *weiPtr = static_cast<PrecisionType *>(this->weightsImageForReconstruction->data);
    for (int v = 0; v < this->correctionImage->nvox; v++)
    {
      corPtr[v] = 0.0;
      weiPtr[v] = 0.0;
    }

    PrecisionType totalSSD = 0.;

    // Iterate over all dynamic iamges
    for (size_t n = 0; n < this->numberOfDynamicImages; ++n)
    {
      // Get the transformation from the current surrogate signal 
      auto curTrafo = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[n] );

      // Get access to the image in acquisition space.
      nifti_image* curDynamicImg = this->dynamicImages[n];

      // Perform the 'forward projection'
      // The forward projection simulates the image acquisition and gets us to the acquisition space
      //
      // Warp the current reconstructed image, then project it into the acquisition (dynamic) space
      nifti_image* curMinSizeWarpedImageInFullImgSpace = this->imageAcquisition->AllocateMinimumSizeImgInFullImgSpace( this->reconstructedImage, curDynamicImg, n );
      curTrafo->TransformImage( this->reconstructedImage, curMinSizeWarpedImageInFullImgSpace );
      nifti_image* simulatedDynamicImage = this->imageAcquisition->SimulateImageAcquisition( curMinSizeWarpedImageInFullImgSpace, curDynamicImg, n );

      // Allocate and calculate the differnce image and SSD between the simulated acquired and acquired iamge
      // The difference image is only used locally, so no need to make it a member of this class.
      nifti_image* differenceImage = nifti_copy_nim_info( curDynamicImg );
      differenceImage->scl_slope = 1.f;
      differenceImage->scl_inter = 0.f;
      differenceImage->nbyper = sizeof( PrecisionType );
      if (sizeof( PrecisionType ) == sizeof( float ))
        differenceImage->datatype = NIFTI_TYPE_FLOAT32;
      else
        differenceImage->datatype = NIFTI_TYPE_FLOAT64;
      
      differenceImage->data = (void*)calloc( differenceImage->nvox, differenceImage->nbyper );

      // Pointers for the SSD calculation
      PrecisionType dynSSD = 0.;
      PrecisionType* dynPtr = static_cast<PrecisionType*>(curDynamicImg->data);
      PrecisionType* simDynPtr = static_cast<PrecisionType*>(simulatedDynamicImage->data);
      PrecisionType* diffPtr = static_cast<PrecisionType*>(differenceImage->data);
      paddingVal = static_cast<PrecisionType>(curTrafo->GetPaddingValue());

      for (size_t vox = 0; vox < differenceImage->nvox; ++vox)
      {
        // Only calculate the difference if both pointer point to values other than nan.
        if (dynPtr[vox] == dynPtr[vox] && simDynPtr[vox] == simDynPtr[vox])
        {
          diffPtr[vox] = dynPtr[vox] - simDynPtr[vox];
          dynSSD += diffPtr[vox] * diffPtr[vox];
        }
        else
        {
          diffPtr[vox] = paddingVal;
        }
      }

      // Add the SSD value evaluated for this dynamic image to the total SSD
      totalSSD += dynSSD;

      // Now bring the difference image into the full image spcae (back-projection)

      // Use the acquisition to calcualte the image after the adjoint operation
      this->imageAcquisition->CalculateAdjoint( this->reconstructionGeometryImage, differenceImage, n );

      // Then transform it using the push interpolation (which accummulates the result in the destination image)
      curTrafo->TransformImageAdjoint( this->imageAcquisition->GetImageAfterAdjoint(), 
        this->imageAcquisition->GetWeightsImageAfterAdjoint(), 
        this->correctionImage, 
        this->weightsImageForReconstruction );

      // Clean up
      nifti_image_free( differenceImage );

      // It is possible that when no image acquisition is simulated the images current simulated dynamic 
      // image and the current minimally sized warped image in full image space are identical. Only need to 
      // free the memory once. 
      if (curMinSizeWarpedImageInFullImgSpace != simulatedDynamicImage)
      {
        nifti_image_free( curMinSizeWarpedImageInFullImgSpace );
        nifti_image_free( simulatedDynamicImage );
      }
      else
      {
        nifti_image_free( curMinSizeWarpedImageInFullImgSpace );
      }
    } // iterate over dynamic iamges

    // normaliset the correction image by weights image
    for (size_t v = 0; v < this->correctionImage->nvox; ++v)
    {
      corPtr[v] = corPtr[v] / weiPtr[v];
    }

    // Check that the reconstruction improved the results, stop itearting, if not
    if ((totalSSD < prevSSD) || (prevSSD == -1))
    {
      // Store the previous reconstruction result
      memcpy( prevReconImage->data, this->reconstructedImage->data, prevReconImage->nvox * prevReconImage->nbyper );
      prevSSD = totalSSD;

      // Now add the correction image to the current reconstruction image. 
      PrecisionType* recPtr = static_cast<PrecisionType*>(this->reconstructedImage->data);
      for (size_t v = 0; v < this->correctionImage->nvox; ++v)
      {
        recPtr[v] = recPtr[v] + corPtr[v];
      }
    }
    else
    {
      break;
    }
  } // reconstruction iterations

  // revert recon (floating image) to previous recon  
  /// ToDoL Figure out exactly why this is.
  memcpy( this->reconstructedImage->data, prevReconImage->data, prevReconImage->nvox * prevReconImage->nbyper );

  //replace voxels in recon image with weight of 0 with padding value
  PrecisionType *weiPtr = static_cast<PrecisionType *>(this->weightsImageForReconstruction->data);
  PrecisionType *reconPtr = static_cast<PrecisionType *>(this->reconstructedImage->data);

  for (size_t v = 0; v < this->weightsImageForReconstruction->nvox; ++v)
  {
    // NOTE: Could possibply be changed to  weiPtr[v] < 1. If a small number of weights were 
    // accumulated in a voxel, such values should not be trusted.
    if (weiPtr[v] == 0)
      reconPtr[v] = paddingVal;
  }

  //clear previous recon image
  nifti_image_free( prevReconImage );
  return;
}




//-------------------------------------
// MoCoReconSuperResolutionIBP::Update
//-------------------------------------
void MoCoReconSuperResolutionIBP::AllocateCorrectionImage()
{
  // Check that the image that defines the geometry exists
  if (nullptr == this->reconstructionGeometryImage)
  {
    supremo_print_error( "Geometry image has to be set before allocating correction image." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  this->ClearCorrectionImage();
  
  // Assure that the correction image corresponds to the precision type
  this->correctionImage = nifti_copy_nim_info( this->reconstructionGeometryImage );
  
  this->correctionImage->scl_slope = 1.f;
  this->correctionImage->scl_inter = 0.f;
  this->correctionImage->nbyper = sizeof( PrecisionType );
  if (sizeof( PrecisionType ) == sizeof( float ))
    this->correctionImage->datatype = NIFTI_TYPE_FLOAT32;
  else
    this->correctionImage->datatype = NIFTI_TYPE_FLOAT64;
  this->correctionImage->data = (void *)calloc( this->correctionImage->nvox, this->correctionImage->nbyper );

  return;
}




//---------------------------------------------------
// MoCoReconSuperResolutionIBP::ClearCorrectionImage
//---------------------------------------------------
void MoCoReconSuperResolutionIBP::ClearCorrectionImage()
{
  // Only clear up if needed
  if (nullptr != this->correctionImage)
  {
    nifti_image_free( this->correctionImage );
    this->correctionImage = nullptr;
  }
}




//-----------------------------------------------------------------
// MoCoReconWeightedAverage::AllocateWeightsImageForReconstruction
//-----------------------------------------------------------------
void MoCoReconSuperResolutionIBP::AllocateWeightsImageForReconstruction()
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
void MoCoReconSuperResolutionIBP::ClearWeightsImageForReconstruction()
{
  // Only clear the weights image if it was not empty before
  if (nullptr != this->weightsImageForReconstruction)
  {
    nifti_image_free( this->weightsImageForReconstruction );
    this->weightsImageForReconstruction = nullptr;
  }
}





bool MoCoReconSuperResolutionIBP::CheckReconstrutedImageCanBeUpdated() 
{
  // Check if the reconstructed image exists. If not it clearly cannot be updated 
  if (nullptr == this->reconstructedImage)
  {
    return false;
  }

  // Check if the geometry of the reconstructed image and the geometry image is still the same
  if (this->reconstructedImage->nvox != this->reconstructionGeometryImage->nvox) return false;
  if (this->reconstructedImage->nx != this->reconstructionGeometryImage->nx) return false;
  if (this->reconstructedImage->ny != this->reconstructionGeometryImage->ny) return false;
  if (this->reconstructedImage->nz != this->reconstructionGeometryImage->nz) return false;
  if (this->reconstructedImage->dx != this->reconstructionGeometryImage->dx) return false;
  if (this->reconstructedImage->dy != this->reconstructionGeometryImage->dy) return false;
  if (this->reconstructedImage->dz != this->reconstructionGeometryImage->dz) return false;
  if (this->reconstructedImage->sform_code != this->reconstructionGeometryImage->sform_code) return false;
  if (this->reconstructedImage->qform_code != this->reconstructionGeometryImage->qform_code) return false;
  
  // Only need to check sform-matrices if sformCode !=0
  if (this->reconstructedImage->sform_code != NIFTI_XFORM_UNKNOWN)
  {
    for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j) 
      {
        if (this->reconstructedImage->sto_ijk.m[i][j] != this->reconstructionGeometryImage->sto_ijk.m[i][j]) return false;
        if (this->reconstructedImage->sto_xyz.m[i][j] != this->reconstructionGeometryImage->sto_xyz.m[i][j]) return false;
      }
    }
  }
  // Only need to check qform-matrices if qformCode !=0
  if (this->reconstructedImage->qform_code != NIFTI_XFORM_UNKNOWN)
  {
    for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
      {
        if (this->reconstructedImage->qto_ijk.m[i][j] != this->reconstructionGeometryImage->qto_ijk.m[i][j]) return false;
        if (this->reconstructedImage->qto_xyz.m[i][j] != this->reconstructionGeometryImage->qto_xyz.m[i][j]) return false;
      }
    }
  }

  return true;
}



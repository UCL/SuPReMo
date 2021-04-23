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
#include "MoCoRecon.h"
#include "Supremo.h"




//----------------------
// MoCoRecon::MoCoRecon
//----------------------
MoCoRecon::MoCoRecon() :
  correspondenceModel( nullptr ),
  imageAcquisition( nullptr ),
  reconstructedImage( nullptr ),
  reconstructionGeometryImage( nullptr ),
  numberOfDynamicImages( 0 )
{}




//-----------------------
// MoCoRecon::~MoCoRecon
//-----------------------
MoCoRecon::~MoCoRecon()
{
  this->ClearReconstructedImage();
  this->imageAcquisition = nullptr;
  this->correspondenceModel = nullptr;
}




//-----------------------------------
// MoCoRecon::SetCorrespondenceModel
//-----------------------------------
void MoCoRecon::SetCorrespondenceModel( const std::shared_ptr<CorrespondenceModel>& correspondenceModelIn )
{
  this->correspondenceModel = correspondenceModelIn;
}




//--------------------------------
// MoCoRecon::SetImageAcquisition
//--------------------------------
void MoCoRecon::SetImageAcquisition( const std::shared_ptr<ImageAcquisition>& imageAcquisitionIn )
{
  this->imageAcquisition = imageAcquisitionIn;
}




//--------------------------------
// MoCoRecon::SetSurrogateSignals
//--------------------------------
void MoCoRecon::SetSurrogateSignals( const SurrogateSignalType & surrSignalsIn )
{
  this->surrogateSignals = surrSignalsIn;
}




//-----------------------------
// MoCoRecon::SetDynamicImages
//-----------------------------
void MoCoRecon::SetDynamicImages( const std::vector<nifti_image*>& dynamicImagesIn )
{
  this->dynamicImages = dynamicImagesIn;
  this->numberOfDynamicImages = this->dynamicImages.size();
}




//-------------------------------------------
// MoCoRecon::SetReconstructionGeometryImage
//-------------------------------------------
void MoCoRecon::SetReconstructionGeometryImage( nifti_image * reconstructionGeometryImageIn )
{
  this->reconstructionGeometryImage = reconstructionGeometryImageIn;
}




//----------------------------------
// MoCoRecon::GetReconstructedImage
//----------------------------------
nifti_image* MoCoRecon::GetReconstructedImage()
{
  return this->reconstructedImage;
}




//--------------------------------------------------
// MoCoRecon::CopyReconstructedImageContentsToImage
//--------------------------------------------------
void MoCoRecon::CopyReconstructedImageContentsToImage( nifti_image * destinationImage )
{
  // Check that reconstructed image and geometry image exist
  if (nullptr == destinationImage || nullptr == this->reconstructedImage)
  {
    supremo_print_error( "Reconstructed image and destination image must not be null" );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Check minimum requirements that the images are equivalent
  if (this->reconstructedImage->nvox     != destinationImage->nvox     ||
      this->reconstructedImage->datatype != destinationImage->datatype ||
      this->reconstructedImage->nbyper   != destinationImage->nbyper)
  {
    supremo_print_error( "Reconstructed image and geometry image must be of same size" );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Copy the data into the destination image
  memcpy( destinationImage->data, this->reconstructedImage->data, this->reconstructedImage->nvox * this->reconstructedImage->nbyper );
}




//---------------------------------------
// MoCoRecon::AllocateReconstructedImage
//---------------------------------------
void MoCoRecon::AllocateReconstructedImage()
{
  // Check that the geometry of the reconstructed image was defined
  if (nullptr == this->reconstructionGeometryImage)
  {
    supremo_print_error( "Reconstruction geometry has to be set before allocating reconstructed image" );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Clear the reconstructed image if it existed before
  this->ClearReconstructedImage();

  // Copy the header information of the geometry-defining image into the reconstructed image
  this->reconstructedImage = nifti_copy_nim_info( this->reconstructionGeometryImage );
  this->reconstructedImage->nbyper = sizeof( PrecisionType );
  
  if (sizeof( PrecisionType ) == sizeof( float ))
    this->reconstructedImage->datatype = NIFTI_TYPE_FLOAT32;
  else
    this->reconstructedImage->datatype = NIFTI_TYPE_FLOAT64;
  
  this->reconstructedImage->data = (void *)calloc( this->reconstructedImage->nvox, this->reconstructedImage->nbyper );
}




//------------------------------------
// MoCoRecon::ClearReconstructedImage
//------------------------------------
void MoCoRecon::ClearReconstructedImage()
{
  if (nullptr != this->reconstructedImage)
  {
    nifti_image_free( this->reconstructedImage );
    this->reconstructedImage = nullptr;
  }
}



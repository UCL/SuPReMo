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




#include "ImageAcquisition.h"
#include "Supremo.h"




//------------------------------------
// ImageAcquisition::ImageAcquisition
//------------------------------------
ImageAcquisition::ImageAcquisition()
  : imageAfterAdjoint(nullptr),
  weightsImageAfterAdjoint(nullptr),
  curImageInAcquisitionSpace(nullptr),
  curImageInFullImgSpace(nullptr)
{}




//-------------------------------------
// ImageAcquisition::~ImageAcquisition
//-------------------------------------
ImageAcquisition::~ImageAcquisition()
{
  // General clean up
  this->ClearImageAfterAdjoint();
  this->ClearWeightsImageAfterAdjoint();
}




//----------------------------------------
// ImageAcquisition::GetImageAfterAdjoint
//----------------------------------------
nifti_image* ImageAcquisition::GetImageAfterAdjoint()
{
  return this->imageAfterAdjoint;
}




//-----------------------------------------------
// ImageAcquisition::GetWeightsImageAfterAdjoint
//-----------------------------------------------
nifti_image* ImageAcquisition::GetWeightsImageAfterAdjoint()
{
  return this->weightsImageAfterAdjoint;
}




//----------------------------------------------------
// ImageAcquisition::AllocateWeightsImageAfterAdjoint
//----------------------------------------------------
void ImageAcquisition::AllocateWeightsImageAfterAdjoint()
{
  //dynamic after adjoint image must already have been allocated
  if (this->imageAfterAdjoint == nullptr)
  {
    supremo_print_error( "Image after adjoint has to be allocated before allocating corresponding weights image." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }
  
  this->ClearWeightsImageAfterAdjoint();
  
  this->weightsImageAfterAdjoint = nifti_copy_nim_info( this->imageAfterAdjoint );
  this->weightsImageAfterAdjoint->scl_slope = 1.f;
  this->weightsImageAfterAdjoint->scl_inter = 0.f;
  this->weightsImageAfterAdjoint->nbyper = sizeof( PrecisionType );
  if (sizeof( PrecisionType ) == sizeof( float ))
    this->weightsImageAfterAdjoint->datatype = NIFTI_TYPE_FLOAT32;
  else
    this->weightsImageAfterAdjoint->datatype = NIFTI_TYPE_FLOAT64;
  this->weightsImageAfterAdjoint->data = (void *)calloc( this->weightsImageAfterAdjoint->nvox, this->weightsImageAfterAdjoint->nbyper );
}




//------------------------------------------
// ImageAcquisition::ClearImageAfterAdjoint
//------------------------------------------
void ImageAcquisition::ClearImageAfterAdjoint()
{
  if (this->imageAfterAdjoint != nullptr)
  {
    nifti_image_free( this->imageAfterAdjoint );
    this->imageAfterAdjoint = nullptr;
  }
}




//-------------------------------------------------
// ImageAcquisition::ClearWeightsImageAfterAdjoint
//-------------------------------------------------
void ImageAcquisition::ClearWeightsImageAfterAdjoint()
{
  if (nullptr != this->weightsImageAfterAdjoint)
  {
    nifti_image_free( this->weightsImageAfterAdjoint );
    this->weightsImageAfterAdjoint = nullptr;
  }
}

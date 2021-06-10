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




#include <iostream>
#include "_reg_tools.h"



 //----------------------------
 // ImagePyramid::ImagePyramid
 //----------------------------
template<class VoxelType>
ImagePyramid<VoxelType>::ImagePyramid() :
  numberOfLevels( 0 ),
  numberOfLevelsToPerform( 0 ),
  imageLevels( nullptr )
{
#ifdef _DEBUG
  std::cout << "Called ImagePyramid<VoxelType>::ImagePyramid()" << std::endl;
#endif
}



//---------------------------------------------------
// ImagePyramid::ImagePyramid( const ImagePyramid& )
//---------------------------------------------------
template<class VoxelType>
ImagePyramid<VoxelType>::ImagePyramid( const ImagePyramid& rhsImagePyramid )
{
  // Copy basic information 
  this->numberOfLevels = rhsImagePyramid->numberOfLevels;
  this->numberOfLevelsToPerform = rhsImagePyramid->numberOfLevelsToPerform;

  // Copy image content of levels

  this->imageLevels = (nifti_image **)malloc( this->numberOfLevelsToPerform * sizeof( nifti_image * ) );

  for (unsigned int iLevel = 0; iLevel < this->numberOfLevelsToPerform; ++iLevel)
  {
    // Copy the nifti-structure if it is not a null pointer
    if (nullptr != rhsImagePyramid->imageLevels[iLevel])
    {
      this->imageLevels[iLevel] = nifti_copy_nim_info( rhsImagePyramid->imageLevels[iLevel] );

      // If the data exists, allocate new memory and copy contents over from the given pyramid 
      if (nullptr != rhsImagePyramid->imageLevels[iLevel]->data)
      {
        this->imageLevels[iLevel]->data = malloc( this->imageLevels[iLevel]->nvox * this->imageLevels[iLevel]->nbyper );
        memcpy( this->imageLevels[iLevel]->data, rhsImagePyramid->imageLevels[iLevel]->data, this->imageLevels[iLevel]->nvox * this->imageLevels[iLevel]->nbyper );
      }
      else 
      {
        this->imageLevels[iLevel]->data = nullptr;
      }
    }
    else
    {
      this->imageLevels[iLevel] = nullptr;
    }
  }
}




//--------------------------
 // ImagePyramid::operator=
 //-------------------------
template<class VoxelType>
ImagePyramid<VoxelType>& ImagePyramid<VoxelType>::operator=( const ImagePyramid<VoxelType> & rhsImagePyramid)
{
  
  // Nothing to copy if the right-hand side is the same as the current object
  if (this == &rhsImagePyramid)
  {
    return *this;
  }

  // Clean up the image in every level before proceeding
  for (unsigned int iLevel = 0; iLevel < this->numberOfLevelsToPerform; ++iLevel)
  {
    if (nullptr != this->imageLevels[iLevel])
    {
      nifti_image_free( this->imageLevels[iLevel] );
      this->imageLevels[iLevel] = nullptr;
    }
  }

  // Clean up the pointer holding the images
  if (this->imageLevels != nullptr) {
    free( this->imageLevels );
    this->imageLevels = nullptr;
  }

  // Copy the level information
  this->numberOfLevels = rhsImagePyramid->numberOfLevels;
  this->numberOfLevelsToPerform = rhsImagePyramid->numberOfLevelsToPerform;

  // And the images itself
  // Allocate memory for the image pointers for all necessary levels
  if (rhsImagePyramid->imageLevels == nullptr)
  {
    this->imageLevels = nullptr;
  }
  else
  {
    this->imageLevels = (nifti_image **)malloc( this->numberOfLevelsToPerform * sizeof( nifti_image * ) );

    for (unsigned int iLevel = 0; iLevel < this->numberOfLevelsToPerform; ++iLevel)
    {
      // copy header information
      this->imageLevels[iLevel] = nifti_copy_nim_info( rhsImagePyramid->GetLevel( iLevel ) );

      // and image contents 
      this->imageLevels[iLevel]->data = malloc( this->imageLevels[iLevel]->nvox * this->imageLevels[iLevel]->nbyper );
      memcpy( this->imageLevels[iLevel]->data, rhsImagePyramid->GetLevel( iLevel )->data, this->imageLevels[iLevel]->nvox * this->imageLevels[iLevel]->nbyper );
    }
  }

  return *this;
}




//------------------------------
 // ImagePyramid::~ImagePyramid
 //-----------------------------
template<class VoxelType>
ImagePyramid<VoxelType>::~ImagePyramid()
{
  // Clean up every level, one at a time
  for ( unsigned int level = 0; level < this->numberOfLevelsToPerform; ++level )
  {
    if ( this->imageLevels[level] != nullptr )
    {
      nifti_image_free( this->imageLevels[level] );
    }
  }

  // Clean up the allocated memory for the pointers to the levels.
  if ( this->imageLevels != nullptr ){
    free( this->imageLevels );
    this->imageLevels = nullptr;
  }

#ifdef _DEBUG
  std::cout << "Called ImagePyramid<VoxelType>::~ImagePyramid()" << std::endl;
#endif
}



//-------------------------------
 // ImagePyramid::GenerateLevels
 //------------------------------
template<class VoxelType>
void ImagePyramid<VoxelType>::GenerateLevels(const nifti_image * const imageIn, 
                                              unsigned int numberOfLevelsIn, 
                                              unsigned int numberOfLevelsToPerformIn, 
                                              unsigned int minNumberOfVoxelsPerDim, 
                                              float( *minVoxDim )[3] )
{
#ifdef _DEBUG
  std::cout << "Called ImagePyramid<VoxelType>::GenerateLevels" << std::endl;
#endif
  
  // Save the input parameters to the member variables
  this->numberOfLevels = (numberOfLevelsToPerformIn > numberOfLevelsIn) ? numberOfLevelsToPerformIn : numberOfLevelsIn;
  this->numberOfLevelsToPerform = numberOfLevelsToPerformIn;
  
  // Allocate memory for the image pointers for all necessary levels
  this->imageLevels = (nifti_image **) malloc( this->numberOfLevelsToPerform * sizeof( nifti_image * ) );

  // First copy over the image information
  // The highest resolution level will be the last that is being processed, i.e. [numberOfLevelsToPerform-1]
  this->imageLevels[this->numberOfLevelsToPerform - 1] = nifti_copy_nim_info( imageIn );

  // Allocate memory for the image data
  this->imageLevels[this->numberOfLevelsToPerform - 1]->data
    = (void*) malloc( this->imageLevels[this->numberOfLevelsToPerform - 1]->nvox *
                      this->imageLevels[this->numberOfLevelsToPerform - 1]->nbyper );

  // Copy over the image contents 
  memcpy( this->imageLevels[this->numberOfLevelsToPerform - 1]->data, imageIn->data,
          this->imageLevels[this->numberOfLevelsToPerform - 1]->nvox * this->imageLevels[this->numberOfLevelsToPerform - 1]->nbyper );

  // Convert to required data format and remove scl-info
  reg_tools_changeDatatype<VoxelType>( this->imageLevels[this->numberOfLevelsToPerform - 1] );
  
  /// \todo Need to find a solution if slope is zero. If this is the case, the input data will be blacked out. 
  // This correction was performed in reg-resp on reading the image
  if (this->imageLevels[this->numberOfLevelsToPerform - 1]->scl_slope == 0)
    this->imageLevels[this->numberOfLevelsToPerform - 1]->scl_slope = 1.0f;

  reg_tools_removeSCLInfo( this->imageLevels[this->numberOfLevelsToPerform - 1] );

  // Need to down-sample the highest resolution image, if not all levels are processed, 
  // i.e. difference between numberOfLevelsToPerform and numberOfLevels
  for ( unsigned int l = this->numberOfLevelsToPerform; l < this->numberOfLevels; l++ )
  {
    bool downsampleAxis[8] = { false, true, true, true, false, false, false, false };
    if ( static_cast<unsigned int>( this->imageLevels[this->numberOfLevelsToPerform - 1]->nx ) <= minNumberOfVoxelsPerDim ) downsampleAxis[1] = false;
    if ( static_cast<unsigned int>( this->imageLevels[this->numberOfLevelsToPerform - 1]->ny ) <= minNumberOfVoxelsPerDim ) downsampleAxis[2] = false;
    if ( static_cast<unsigned int>( this->imageLevels[this->numberOfLevelsToPerform - 1]->nz ) <= minNumberOfVoxelsPerDim ) downsampleAxis[3] = false;

    //if minVoxDim provided then only downsample if resampled vox in image <= minVoxDim
    if ( minVoxDim != nullptr )
    {
      if ( this->imageLevels[this->numberOfLevelsToPerform - 1]->dx * 2.f > minVoxDim[this->numberOfLevelsToPerform - 1][0] ) downsampleAxis[1] = false;
      if ( this->imageLevels[this->numberOfLevelsToPerform - 1]->dy * 2.f > minVoxDim[this->numberOfLevelsToPerform - 1][1] ) downsampleAxis[2] = false;
      if ( this->imageLevels[this->numberOfLevelsToPerform - 1]->dz * 2.f > minVoxDim[this->numberOfLevelsToPerform - 1][2] ) downsampleAxis[3] = false;
    }

    reg_downsampleImage<VoxelType>( this->imageLevels[this->numberOfLevelsToPerform - 1], 1, downsampleAxis );
  }

  // Now generate the next layers
  // Image for each subsequent level is allocated and downsampled
  for ( int l = this->numberOfLevelsToPerform - 2; l >= 0; --l )
  {
    // Copy the image information
    this->imageLevels[l] = nifti_copy_nim_info( this->imageLevels[l + 1] );
    
    // Allocate memory and copy image contents
    this->imageLevels[l]->data = (void *) calloc( this->imageLevels[l]->nvox,
                                                  this->imageLevels[l]->nbyper );
    memcpy( this->imageLevels[l]->data, 
            this->imageLevels[l + 1]->data,
            this->imageLevels[l]->nvox * this->imageLevels[l]->nbyper );

    // Downsample the image (if appropriate)
    bool downsampleAxis[8] = { false, true, true, true, false, false, false, false };
    if ( static_cast<unsigned int>( this->imageLevels[l]->nx ) <= minNumberOfVoxelsPerDim ) downsampleAxis[1] = false;
    if ( static_cast<unsigned int>( this->imageLevels[l]->ny ) <= minNumberOfVoxelsPerDim ) downsampleAxis[2] = false;
    if ( static_cast<unsigned int>( this->imageLevels[l]->nz ) <= minNumberOfVoxelsPerDim ) downsampleAxis[3] = false;
    
    // if minVoxDim provided then only downsample if resampled vox in image <= minVoxDim
    if ( minVoxDim != nullptr )
    {
      if ( this->imageLevels[l]->dx * 2 > minVoxDim[l][0] ) downsampleAxis[1] = false;
      if ( this->imageLevels[l]->dy * 2 > minVoxDim[l][1] ) downsampleAxis[2] = false;
      if ( this->imageLevels[l]->dz * 2 > minVoxDim[l][2] ) downsampleAxis[3] = false;
    }

    reg_downsampleImage<VoxelType>( this->imageLevels[l], 1, downsampleAxis );
  }
}




//-------------------------
 // ImagePyramid::GetLevel
 //------------------------
template<class VoxelType>
nifti_image* ImagePyramid<VoxelType>::GetLevel( unsigned int level )
{

  // Check that the requested level is within the current range
  if ( (level < 0) || (level >= this->numberOfLevelsToPerform) )
  {

#ifdef _DEBUG
    std::cout << "WARNING: Requested level does not exist!" << std::endl;
#endif

    return  nullptr;
  }

  return this->imageLevels[level];
}

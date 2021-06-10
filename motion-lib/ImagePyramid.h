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


#pragma once


#include "nifti1.h" 
#include "nifti1_io.h"


/** Class that implements a mulit-resolution image pyramid. 
* This class implements a decimating multi-resolution image pyramid. From one level to the 
* next the image resolution is halved. 
*/
template<class VoxelType>
class ImagePyramid{

public:
  /** Constructor 
   */
  ImagePyramid();
  
  /** Copy constructor
   */
  ImagePyramid( const ImagePyramid& rhsPyramid );
  
  /** Destructor. Deletes the image data that was generated with GenerateLevels() 
   */
  ~ImagePyramid();

  /** Copy operator
   */
  ImagePyramid & operator= ( const ImagePyramid<VoxelType>& rhsPyramid );


  /** Generate the pyramid levels
   * \param imageIn Pointer to the nifti input image
   * \param numberOfLevels Total number of pyramid levels. Will be increased, if numberOfLevelsToPerform is larger than this.
   * \param numberOfLevelsToPerform Total number of levels to perform. If it is larger than numberOfLevels, it will be updated.
   * \param minNumberOfVoxelsPerDim The minimum number of voxels per level below which a dimension will not be subsampled.
   * \param minVoxDim The minimal voxel dimension. If minVoxDim is provided, then the image is only downsampled if 
                      voxel size in the resampled image <= minVoxDim.
   */
  void GenerateLevels( const nifti_image* const imageIn,
                       unsigned int numberOfLevels,
                       unsigned int numberOfLevelsToPerform,
                       unsigned int minNumberOfVoxelsPerDim = 1, 
                       float (* minVoxDim)[3] = nullptr );
  
  /** Get the pointer to a nifti image that represents a level of the decimating Gaussian image pyramid.
   * \param level The requested level. 
   * \return A pointer to the requested pyramid level. Will return NULL if the requested level does not exist.
   */
  inline nifti_image* GetLevel(unsigned int level);




private:
  unsigned int numberOfLevels;           ///< The total number of levels
  unsigned int numberOfLevelsToPerform;  ///< The number of levels on which computations are performed

  nifti_image ** imageLevels;            ///< Pointer to all image levels
};


#include "ImagePyramid.cpp"
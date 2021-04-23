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

#include "nifti1_io.h"

/** Base class defining the image similarity interface
 */
class ImageSimilarity{
public:
  /** Measures and returns the image similarity \f$ S(\mathbf{P}_t, \mathbf{P}_{A_t}) \f$ between two images. 
   */
  virtual double GetSimilarityMeasureValueForImages( nifti_image* img1In,
                                                     nifti_image* img2In ) = 0;
  
  /** Calculates and returns the gradient of the image similarity measure with respect to the voxel intensities.
   */
  virtual void GetSimilarityGradientWRTVoxels( nifti_image* referenceImg,
											                         nifti_image* sourceImg,
                                               nifti_image* similarityGradientWRTVoxelsOutputImage ) = 0;
};

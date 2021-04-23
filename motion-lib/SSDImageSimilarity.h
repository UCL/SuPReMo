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

#include "ImageSimilarity.h"


/** Class to calculate the sum-of-squared-differences (SSD) image similarity value and gradient.
 */
class SSDImageSimilarity:public ImageSimilarity
{
public:
  typedef float PrecisionType;

  /** Constructor.
   */
  SSDImageSimilarity();


  /** Destructor.
   */
  ~SSDImageSimilarity();
  
  /** Measures and returns the image similarity \f$ S(\mathbf{P}_t, \mathbf{P}_{A_t}) \f$ between two images. 
   * 
   */
  double GetSimilarityMeasureValueForImages( nifti_image* source,
                                             nifti_image* reference );
  
  /** Calculates and returns the gradient of the image similarity measure with respect to the voxel intensities.
   *  Unlike SSDImageSimilarity::GetSimilarityMeasureValueForImages this funciton is not symmetric.
   *  Here \f$ \frac{\partial \mathcal{S}_t}{\mathbf{P}_{A_t}}\f$ is calculated, which for SSD is given by 
   *  \f[ \frac{ \partial \mathcal{S} }{ \partial p_{A_t ,y} } = -2\left( p_{A_t ,y} - \partial p_{t ,y} \right) \f] 
   *  calculated at each voxel. 
   *  
   *  \param referenceImg Nifti image pointer to the reference image. In the case of the unified registration framework it is
   *                   the actual dynamic image \f$ \mathbf{P}_{t} \f$.
   *  \param sourceImg Nifti image pointer to the source image. In the case of the unified registration framework it is 
   *                   the simulated dynamic image \f$ \mathbf{P}_{A_t} \f$.
   *  \param similarityGradientWRTVoxelsOutputImage Nifti image of the same size/geometry of the sourceImg and referenceImg
   *                                                that is used to recieve the gradient image. Purposely not a member of this
   *                                                class.
   */
  virtual void GetSimilarityGradientWRTVoxels( nifti_image* referenceImg,
	                                             nifti_image* sourceImg,
                                               nifti_image* similarityGradientWRTVoxelsOutputImage );
private:
  double GetSSDValue( nifti_image* referenceImage, nifti_image* sourceImage );
};



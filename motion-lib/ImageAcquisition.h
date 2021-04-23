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

/** Base class defining the image acquisition interface. 
 */
class ImageAcquisition
{
public:
  typedef float PrecisionType;

  /** Constructor
   */
  ImageAcquisition();
  
  /** Destructor
   */
  virtual ~ImageAcquisition();
   
  /** Simulate the image acquisition. Needs to be implemented by derived class depending on the image acquisition method.
   */
  virtual nifti_image* SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace ) = 0;
  
  /** Calculate the adjoint of the image acquisition. Needs to be implementd by derived class according to the image acquisition method.
   */
  virtual void CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace ) = 0;
  
  /** Allocate an empty image in full image space that has the minimum required size to simulate the image acquisition. 
   *  This allows, for instance, efficient warping of the image in full image space without the need to warp the full 
   *  sized image. This has to be implemented by a derived class. This was usually used as the "warped" image in reg-resp.
   *  \param imgInFullImgSpace The full sized image
   *  \param imgInAcquisitionSpace The image in acquisition space  
   *  \todo Decide how to handle the image deallocation
   */
  virtual nifti_image* AllocateMinimumSizeImgInFullImgSpace( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace ) = 0;
  
  /** Get a pointer to the image after adjoint
   */
  virtual nifti_image* GetImageAfterAdjoint();
  
  /** Get a copy to the weights image after adjoint
   */
  virtual nifti_image* GetWeightsImageAfterAdjoint();

protected:
  
  virtual void AllocateImageAfterAdjoint() = 0;    
  virtual void AllocateWeightsImageAfterAdjoint(); // Can go into base-class since this takes the image after adjoint as a basis for copying the data
  virtual void ClearImageAfterAdjoint();
  virtual void ClearWeightsImageAfterAdjoint();

  nifti_image* imageAfterAdjoint;
  nifti_image* weightsImageAfterAdjoint;
  
  nifti_image* curImageInAcquisitionSpace;
  nifti_image* curImageInFullImgSpace;
};

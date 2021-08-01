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
   
  /** Simulate the image acquisition also known as a forward-projection depending on the imaging modality. Needs to be implemented by derived class 
   *  depending on the image acquisition method. The returned image should be allocated by this function and deleted after use externally(see \c return ). 
   *
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space. 
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space. 
   *  \param dynamicImageTimePoint The dynamic image time point in case time-dependent meta data needs to be used internally.
   *  \return Pointer to the simulated dyanmic image. This image will be allocated by this function and is expected to be deleted after usaage using \c nifti_image_free() . 
   */
  virtual nifti_image* SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace, unsigned int dynamicImageTimePoint ) = 0;
  
  /** Calculate the adjoint of the image acquisition - also known as the back-projection for some imaging modalities. Needs to be implementd by derived 
   *  class according to the image acquisition method. The adjoint image will be allocated internally and kept as a member variable. To access it, use 
   *  \ref GetImageAfterAdjoint. This image will be deleted by this object. 
   * 
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space. 
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space. 
   *  \param dynamicImageTimePoint The dynamic image time point in case time-dependent meta data needs to be used internally.
   */
  virtual void CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint) = 0;
  
  /** Allocate an empty image in full image space that has the minimum required size to simulate the image acquisition. 
   *  This allows, for instance, efficient warping of the image in full image space without the need to warp the full 
   *  sized image. This has to be implemented by a derived class. This was usually used as the "warped" image in reg-resp.
   * 
   *  \param imgInFullImgSpace The full sized image
   *  \param imgInAcquisitionSpace The image in acquisition space  
   *  \param dynamicImageTimePoint Time point of the dynamic image (if used).
   */
  virtual nifti_image* AllocateMinimumSizeImgInFullImgSpace( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint) = 0;
  
  /** Get a pointer to the image after adjoint (after the back-projection)
   *  \return Returns a pointer to the nifti image calculated by the function \ref CalculateAdjoint. 
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
  
  nifti_image* curImageInAcquisitionSpace;  ///< Current image in acquisition space
  nifti_image* curImageInFullImgSpace;      ///< Current image in full image space
  unsigned int curDynamicImageTimePoint;    ///< Current dynamic image time point
};

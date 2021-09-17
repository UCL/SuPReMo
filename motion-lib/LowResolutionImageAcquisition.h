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

#include "ImageAcquisition.h"
#include "Supremo.h"

class LowResolutionImageAcquisition : public ImageAcquisition
{
public:
  typedef float PrecisionType;

  /** Default constructor
   */
  LowResolutionImageAcquisition();

  /** Destructor
   */
  virtual ~LowResolutionImageAcquisition();
   
  /** Implementation of the low-resolution image acquisition simulation. Calculates a low-resolution image using the information provided by the 
   *  image in image space. The simulated image will be allocated and returned as a pointer. The image in acuqisition space will be used to determine
   *  the smoothing and downwampling. 
   * 
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space. 
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space. 
   *  \param dynamicImageTimePoint The dynamic image time point in case time-dependent meta data needs to be used internally.
   *  \return Pointer to the simulated dyanmic image. This image will be allocated by this function and is expected to be deleted after usaage using \c nifti_image_free() . 
   */
  nifti_image* SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace, unsigned int dynamicImageTimePoint);

  /** Backproject the low-resoution image into full image space.
   * 
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space. 
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space. 
   *  \param dynamicImageTimePoint Not used in this method.
   */
  void CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint);

  /** Computes the minimum sized image in full image space for efficient warping in the objective function calculations. Only the geometric information of
   *  the images is used, the dynamic image time point is ignored here. 
   * 
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space.
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space.
   *  \param dynamicImageTimePoint Not used in this method.
   */
  nifti_image* AllocateMinimumSizeImgInFullImgSpace( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint);

protected:
  virtual void AllocateImageAfterAdjoint();  
  
  const PrecisionType lowResolutionThreshold;          ///< Ratio between the high and low resolution image. Only if ratio is larger, low resolution acquisition will be simulated
  const PrecisionType roundErrorThreshold;             ///< Value used below which rounding can be savely done without loss of accuracy
  PrecisionType paddingValue;                          ///< The padding value used in simulating the image acquisition
};

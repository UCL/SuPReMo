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

class NoImageAcquisition : public ImageAcquisition
{
public:
  typedef float PrecisionType;

  /** Default constructor
   */
  NoImageAcquisition();

  /** Destructor
   */
  virtual ~NoImageAcquisition();
   
  /** Implementation of the image acquisition simulation. Since no image acqusition is simulated, the image in full space is simply returned. This needs
   *  to be considered when deleting images, check for equality. 
   *  \param imgInFullImgSpace Pointer to nifti image holding the image in full image space. 
   *  \param imgInAcquisitionSpace Pointer to the nifti image in acquisition space. 
   *  \param dynamicImageTimePoint The dynamic image time point in case time-dependent meta data needs to be used internally.
   */
  nifti_image* SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace, unsigned int dynamicImageTimePoint );

  /** To achieve the required behaviour  when not simulating any image acquisition, copy current image in acquisition space into the image after applying the adjoint.
   *  Image after adjoint still needs to be accessed via the corresponding getter function. 
   *  \param imgInFullImgSpace The full sized image
   *  \param imgInAcquisitionSpace The image in acquisition space  
   *  \param dynamicImageTimePoint Time point of the dynamic image (not used here).
   */
  void CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint);

  /** Implementation of allocating the minimum-sized image in full image space.  Here no image acquisition is simulated, hence the allocated image 
   *  has to have the size of the image in acquisition space. 
   * 
   *  \param imgInFullImgSpace The full sized image
   *  \param imgInAcquisitionSpace The image in acquisition space  
   *  \param dynamicImageTimePoint Time point of the dynamic image (not used).
   *  \return An empty image with the size of \c imgInAcquisitionSpace
   */
  nifti_image* AllocateMinimumSizeImgInFullImgSpace( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace, unsigned int dynamicImageTimePoint);

protected:
  virtual void AllocateImageAfterAdjoint();  
};

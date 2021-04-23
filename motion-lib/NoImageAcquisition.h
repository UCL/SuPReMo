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
   
  /** Implementation of the image acquisition simulation. 
   */
  nifti_image* SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace );

  void CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace );

  nifti_image* AllocateMinimumSizeImgInFullImgSpace( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace );

protected:
  virtual void AllocateImageAfterAdjoint();  
};

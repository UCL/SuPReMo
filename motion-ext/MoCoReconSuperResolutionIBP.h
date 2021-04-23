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

#include "MoCoRecon.h"


/** This function performs up to maxIter iterations of a motion compensated iterative backprojection super-resolution reconstruction. 
 *  It will attempt to minimise the SSD between the simulated acquired image data and the real acquired image data.
 */
class MoCoReconSuperResolutionIBP : public MoCoRecon
{
public:
  // Constructor
  MoCoReconSuperResolutionIBP( unsigned int maxIterationNumber, bool updateReconstruction );

  // Destructor
  virtual ~MoCoReconSuperResolutionIBP();

  /// Run the motion-compensated image reconstruction
  virtual void Update();

private:

  void AllocateWeightsImageForReconstruction();
  void ClearWeightsImageForReconstruction();

  void AllocateCorrectionImage();
  void ClearCorrectionImage();

  bool CheckReconstrutedImageCanBeUpdated();

  unsigned int maxIterationNumber;             ///< The number of iterations performed during reconstruction. Set during construction. 
  bool updateReconstruction;                   ///< Update the reconstruction instead of restarting a new one from scratch with every update.
  nifti_image* correctionImage;                ///< The correction image that will be used to update the reconstructed image.
  nifti_image* weightsImageForReconstruction;  ///< The weights image used during the reconstruciton process reconstruction.
};

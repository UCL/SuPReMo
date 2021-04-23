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


/** Class implementing the weighted averaging motion-compensated image reconstruction.
 */
class MoCoReconWeightedAverage : public MoCoRecon
{
public:
  /** Constructor
   */
  MoCoReconWeightedAverage();

  /** Destructor
   */
  virtual ~MoCoReconWeightedAverage();

  /** Calculate the motion-compensated image reconstruction
   */
  void Update();


private:
  /** Allocate the memory holding the weights image
   */
  void AllocateWeightsImageForReconstruction();
  
  /** Clean up the weights image 
   */
  void ClearWeightsImageForReconstruction();

  nifti_image* weightsImageForReconstruction;   ///< The weights image
};

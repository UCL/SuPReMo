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

// Can only use sprintf_s, strcpy_s on Windows
#ifndef _WIN32
#define sprintf_s sprintf
#define strcpy_s strcpy
#endif

/**
 * Motion compensated image reconstruction methods
 */
enum t_motionCompensatedReconstruction
{
  NO_RECONSTRUCTION = 0,     ///< no motion compensated reconstruction - use provided static image
  WEIGHTED_AVERAGING,        ///< take weighted average of deformed dynamic images
  SUPER_RESOLUTION_RESTART,  ///< super resolution reconstruction using iterative back-projection method, restart the recon
  SUPER_RESOLUTION_UPDATE    ///< super resolution reconstruction using iterative back-projection method, update existing recon
};


/**
 * Dynamic data types
 */
enum t_dynamicData
{
  SAME_RES_AS_STATIC = 0, ///< full or partial images (slices/slabs/volumes) with same res as static volume, i.e. no image acquisition simulated
  LOWER_RES_THAN_STATIC   ///< full or partial images (slices/slabs/volumes) with lower res than static volume, i.e. simulate acquisition of lower res data using Gaussian PSF
};


/**
 * Transformation types
 */
enum t_transformationType
{
  STANDARD_B_SPLINE = 0,
  SLIDING_B_SPLINE
};



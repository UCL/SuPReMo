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

#include <math.h>

// Functions dedicated to replace reg_math macros
template<typename T>
inline int nmm_floorInt( const T& val )
{
  // This is how it should be 
  return static_cast<int>(floor( val ));
}



// #define nmm_floorInt( a ) ((a) > 0 ? (int)(a) : (int)((a)-1))
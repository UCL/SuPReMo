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

#include "Optimiser.h"


/** Class implementing the conjugate gradient optimiser according to 
 *  nifty-reg-resp and Numerical Recipes in C++. 
 */
class ConjugateGradientOptimiser :
  public Optimiser
{
public:
  
  /** Constructor
   */
  ConjugateGradientOptimiser();


  /** Destructor
   */
  ~ConjugateGradientOptimiser();
  

  /** Initiate the optimisation. The iteration numbers as well as the stored objective function values will be reset. 
   *  \param startingPoint Point where to start the optimisation from. Can be null-pointer to start form zero.
   */
  void Optimise( PrecisionType* startingPoint );

  /** Initialise the optimiser. Initialises all parameters of the Optimiser base class and conjugate gradient optimiser.
   */
  virtual void Initialise();


  /** Set if the gradient sequence should be continued even if the objective function changed.
   *  Reg-resp implemented this to be set to true, but it may be beneficial to set this to false and start with a 
   *  fresh gradient sequence if the motion-compensated image and thus the objective function changed with each 
   *  switch iteration.
   *  \param continueOptimisationIn Define the continuation of gradient calculation. Set to false for a fresh start.
   */
  void SetContinueOptimisation( bool continueOptimisationIn ) { this->continueOptimisation = continueOptimisationIn; };

private:
  /** Calculate the conjugate gradient update
   */
  void CalculateGradientUpdate();
  
  /** Function that save the current gradient to nifti images if the output file name is not empty.
   *  \note This will trigger recalculation of the gradient and will affect performance. 
   */
  void SaveIntermediateGradientImage();

  PrecisionType* array1;                       ///< The array g in numerical recipes
  PrecisionType* array2;                       ///< The array h in numerical recipes

  bool continuationOfOptimisationPossible;     ///< Indicating if array1 and array2 were initialised properly and can be used to continue with the optimisation
  bool continueOptimisation;                   ///< Indicates if the sequence of gradients should be continued. Reg-resp implemented this to be set to true, but it may be beneficial to set this to false and start with a fresh gradient sequence if the motion-compensated image and thus the objective function changed with each switch iteration.
};


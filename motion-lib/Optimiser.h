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

#include <memory>
#include <iostream>

// forward declarations
class ObjectiveFunction;


/** Base class that defines the basic optimiser functionality. Has to be derived and implemented into a working optimiser class.
 */
class Optimiser
{
public:
  typedef float PrecisionType;
  
  /** Constructor
   */
  Optimiser();

  /** Destructor
   */
  virtual ~Optimiser();
  
  /** Set the objective function which is to be optimised
   *  \param objectiveFunctionIn Shared pointer to the objective function
   */
  void SetObjectiveFunction( std::shared_ptr<ObjectiveFunction>& objectiveFunctionIn );
  
  /** Initiate the optimisation
   *  \param startingPoint Point where to start the optimisation from. Can be null-pointer to start form zero. 
   */
  virtual void Optimise( PrecisionType* startingPoint ) = 0;
  
  /** Initialise the optimiser. Needs to be implemented by inheriting class.
   */
  virtual void Initialise() = 0;

  /** Return the best objecitve function value 
   */
  PrecisionType GetBestObjectiveFunctionValue();
  
  /** Return the initial objecitve function value
  */
  PrecisionType GetInitialObjectiveFunctionValue();

  /** Return the pointer to the point which resulted in the best functcion value.
   */
  PrecisionType* GetBestPoint();

  /** Set the current iteration number to zero. 
   */
  void ResetCurrentIterationNumber();

  /** Set the number of iterations 
   *  \param maxIterations The maximum number of iterations performed. 
   */
  inline void SetMaxIterations( unsigned int maxIterations ) 
  {
    this->maxNumberOfIterations = maxIterations;
  }

  /** Set a previously evaluated objective function value such that the first value calculation
   *  can be skipped during the optimisation. Note this value has to be the same as for the starting
   *  point of the optimisation. No additional checks are being performed internally. 
   */
  void SetPreviouslyEvaluatedObjectiveFunctionValue( PrecisionType objectiveFunctionValue );


  void SetOutputIntermediateGradientFolder( const std::string & outputFolderIn, const std::string & outPrefix ) 
  {
    this->outputIntermediateGradientFolder = outputFolderIn;
    this->outputIntermediateGradientPrefix = outPrefix;
  }


protected:
  
  /** Move the current point so that it moves by step size from the best point into the direction of the gradient.
   *  \f[  \mathbf{p}_\mathrm{cur} = \mathbf{p}_\mathrm{best} + s \cdot \nabla\mathbf{f} \f]
   */
  void inline MoveCurrentPointAlongGradient( PrecisionType stepSize );
  
  /** Perform a line search. Walk along the gradient direction with increasing step sizes. If no more improvements are found, 
   *  continue with halved step sizes until the minimum size was reached, or maximum number of maxNumberOfLineIterations was exceeded.
   *
   *  \param maxLength   The maximum size to which the line-search steps are limited. 
   *  \param smallLength The length below which no more line-search steps are performed.
   *  \param startLength The initial step-size. Will be updated with the total stpe size which can be utilised for consecutive 
   *                     line-searches (i.e. with updated gradients).
   */
  void PerformLineSearch( PrecisionType maxLength, PrecisionType smallLength, PrecisionType& startLength );
  
  PrecisionType currentObjectiveFunctionValue;           ///< The current objective function value corresponding to the current iteration
  PrecisionType bestObjectiveFunctionValue;              ///< The best objective function value achieved
  PrecisionType initialObjectiveFunctionValue;           ///< The objective function value at the starting point of the optimisation
  unsigned int currentIterationNumber;                   ///< Tracker of the current iteration number
  unsigned int maxNumberOfIterations;                    ///< The maximum number of allowed iterations 
  const unsigned int maxNumberOfLineIterations;          ///< The maximum number of iterations used for the line-serach along the gradient
  PrecisionType maxLineSearchStepSize;                   ///< The maximum step size used for the line-search. Has to be initialised in derived class.
  PrecisionType minLineSearchStepSize;                   ///< The minimum step size used for the line-search. Has to be initialised in derived class.
  unsigned int numberOfIterationsPerformed;              ///< The total number of iterations performed
  unsigned int numberOfDOFs;                             ///< The number of degrees of freedom.
  std::shared_ptr<ObjectiveFunction> objectiveFunction;  ///< Shared pointer to the objective function.
  PrecisionType* gradient;                               ///< The array xi in numerical recipes.
  PrecisionType* currentPoint;                           ///< The current point, p in numerical recipes.
  PrecisionType* bestPoint;                              ///< The best point achieved, p in numerical recipes.
  PrecisionType preEvaluatedObjectiveFunctionValue;      ///< A pre-calculated objective function value that will be used during the optimisation if usePreEvaluatedEvaluatedObjectiveFunctionValue is set to true.
  bool usePreEvaluatedEvaluatedObjectiveFunctionValue;   ///< Indicator if the first objective function value calculation can be skipped.
  std::string outputIntermediateGradientFolder;          ///< Folder to which the intermediated objective function gradients will be saved.
  std::string outputIntermediateGradientPrefix;          ///< String holding additional information to differentiate the output from different levels
};

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





//----------
// Includes
//----------
#include "Optimiser.h"
#include "ObjectiveFunction.h"




//-----------------------
// Optimiser::Optimiser
//-----------------------
Optimiser::Optimiser()
  : bestObjectiveFunctionValue( 0 ),
  currentObjectiveFunctionValue( 0 ),
  initialObjectiveFunctionValue( 0 ),
  currentIterationNumber( 0 ),
  maxNumberOfIterations( 0 ),
  maxNumberOfLineIterations( 12 ),
  maxLineSearchStepSize( 0 ),
  numberOfIterationsPerformed( 0 ),
  numberOfDOFs( 0 ),
  objectiveFunction( nullptr ),
  gradient( nullptr ),
  currentPoint( nullptr ),
  bestPoint( nullptr ),
  preEvaluatedObjectiveFunctionValue( 0 ),
  usePreEvaluatedEvaluatedObjectiveFunctionValue( false ),
  outputIntermediateGradientFolder( "" ),
  outputIntermediateGradientPrefix( "" )
{
}




//-----------------------
// Optimiser::~Optimiser
//-----------------------
Optimiser::~Optimiser()
{
  if (this->gradient != nullptr) free( this->gradient );
  if (this->currentPoint != nullptr) free( this->currentPoint );
  if (this->bestPoint != nullptr) free( this->bestPoint );
}




//---------------------------------
// Optimiser::SetObjectiveFunction
//---------------------------------
void Optimiser::SetObjectiveFunction( std::shared_ptr<ObjectiveFunction>& objectiveFunctionIn)
{
	this->objectiveFunction = objectiveFunctionIn;
}




//------------------------------------------
// Optimiser::GetBestObjectiveFunctionValue
//------------------------------------------
Optimiser::PrecisionType Optimiser::GetBestObjectiveFunctionValue()
{
  return this->bestObjectiveFunctionValue;
}





//------------------------------------------
// Optimiser::GetInitialObjectiveFunctionValue
//------------------------------------------
Optimiser::PrecisionType Optimiser::GetInitialObjectiveFunctionValue()
{
  return this->initialObjectiveFunctionValue;
}



//-------------------------
// Optimiser::GetBestPoint
//-------------------------
Optimiser::PrecisionType * Optimiser::GetBestPoint()
{
  return this->bestPoint;
}




//----------------------------------------
// Optimiser::ResetCurrentIterationNumber
//----------------------------------------
void Optimiser::ResetCurrentIterationNumber()
{
  this->currentIterationNumber = 0;
  return;
}




//---------------------------------------------------------
// Optimiser::SetPreviouslyEvaluatedObjectiveFunctionValue
//---------------------------------------------------------
void Optimiser::SetPreviouslyEvaluatedObjectiveFunctionValue( PrecisionType objectiveFunctionValueIn )
{
  this->preEvaluatedObjectiveFunctionValue = objectiveFunctionValueIn;
  this->usePreEvaluatedEvaluatedObjectiveFunctionValue = true;
}




//------------------------------------------
// Optimiser::MoveCurrentPointAlongGradient
//------------------------------------------
void Optimiser::MoveCurrentPointAlongGradient( PrecisionType stepSize ) 
{
  /// \todo Make parallel
  for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
  {
    this->currentPoint[iDOF] = this->bestPoint[iDOF] + stepSize * this->gradient[iDOF];
  }

  return;
}




//------------------------------
// Optimiser::PerformLineSearch
//------------------------------
void Optimiser::PerformLineSearch( PrecisionType maxLength, PrecisionType smallLength, PrecisionType& startLength )
{

  // Ensure that the start length is not above the maximum allowed length
  startLength = startLength > maxLength ? maxLength : startLength;

  unsigned int lineIteration = 0;
  PrecisionType cumulativeLength = 0;
  PrecisionType currentLength = startLength;


  // Keep moving in slightly larger steps along gradient until no improvement or max iterations reached
  while (this->currentIterationNumber < this->maxNumberOfIterations && lineIteration < maxNumberOfLineIterations)
  {
    //move along gradient with step size of "current length"
    this->MoveCurrentPointAlongGradient( -currentLength );

    this->currentObjectiveFunctionValue = this->objectiveFunction->GetValue( this->currentPoint );
    this->currentIterationNumber++;
    lineIteration++;

    // Check if step improved objective function value
    if (this->currentObjectiveFunctionValue > this->bestObjectiveFunctionValue)
    {
#ifndef NDEBUG
      printf( "[Supremo DEBUG] [%i] objective function: %g | Increment %g | ACCEPTED\n",
        (int)this->currentIterationNumber,
        this->currentObjectiveFunctionValue,
        currentLength );

#endif
      // Update the objective function value and best parameter point
      this->bestObjectiveFunctionValue = this->currentObjectiveFunctionValue;
      memcpy( this->bestPoint, this->currentPoint, this->numberOfDOFs * sizeof( PrecisionType ) );

      // Update the total added length
      cumulativeLength += currentLength;

      // Increase the step size
      currentLength *= 1.1f;
      currentLength = (currentLength < maxLength) ? currentLength : maxLength;
    }
    else // objective function value did not improve with this step
    {
#ifndef NDEBUG
      printf( "[Supremo DEBUG] [%i] objective function: %g | Increment %g | REJECTED\n",
        (int)this->currentIterationNumber,
        this->currentObjectiveFunctionValue,
        currentLength );
#endif
      // if objective function value did not improve, break out of loop
      break;
    }
  }

  //now keep halving step size and check if one step along gradient improves objective function
  //repeat until step size < min step size or max iterations reached

  //half step size
  currentLength *= 0.5;

  while (currentLength > smallLength && this->currentIterationNumber < this->maxNumberOfIterations && lineIteration < maxNumberOfLineIterations)
  {
    //move along gradient
    this->MoveCurrentPointAlongGradient( -currentLength );

    this->currentObjectiveFunctionValue = this->objectiveFunction->GetValue( this->currentPoint );
    this->currentIterationNumber++;
    lineIteration++;

    // check for improved objective function value
    if (this->currentObjectiveFunctionValue > this->bestObjectiveFunctionValue)
    {
#ifndef NDEBUG
      printf( "[Supremo DEBUG] [%i] objective function: %g | Increment %g | ACCEPTED\n",
        (int)this->currentIterationNumber,
        this->currentObjectiveFunctionValue,
        currentLength );
#endif
      // Update the objective function value and best parameter point
      this->bestObjectiveFunctionValue = this->currentObjectiveFunctionValue;
      memcpy( this->bestPoint, this->currentPoint, this->numberOfDOFs * sizeof( PrecisionType ) );

      // Update the total added length
      cumulativeLength += currentLength;
    }
    else // otherwise reject the step
    {
#ifndef NDEBUG
      printf( "[NiftyReg DEBUG] [%i] objective function: %g | Increment %g | REJECTED\n",
        (int)this->currentIterationNumber,
        this->currentObjectiveFunctionValue,
        currentLength );
#endif
    }
    // half step size
    currentLength *= 0.5;
  }

  // update the current size for the next iteration
  startLength = cumulativeLength;

  // Set the current parameters to the best parameters
  memcpy( this->currentPoint, this->bestPoint, this->numberOfDOFs * sizeof( PrecisionType ) );
}

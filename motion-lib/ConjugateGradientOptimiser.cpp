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
#include "ConjugateGradientOptimiser.h"
#include "Supremo.h"
#include "_reg_ReadWriteImage.h"
#include <sstream>
#include <iomanip>  



//--------------------------------------------------------
// ConjugateGradientOptimiser::ConjugateGradientOptimiser
//--------------------------------------------------------
ConjugateGradientOptimiser::ConjugateGradientOptimiser()
  : array1( nullptr ),
  array2( nullptr ),
  continuationOfOptimisationPossible( false ),
  continueOptimisation(true)
{}




//---------------------------------------------------------
// ConjugateGradientOptimiser::~ConjugateGradientOptimiser
//---------------------------------------------------------
ConjugateGradientOptimiser::~ConjugateGradientOptimiser()
{
  // Clean up the arrays allocated during initialisation
  if (this->array1 != nullptr) free( this->array1 );
  if (this->array2 != nullptr) free( this->array2 );
}




//--------------------------------------
// ConjugateGradientOptimiser::Optimise
//--------------------------------------
void ConjugateGradientOptimiser::Optimise( PrecisionType* startingPoint )
{
  // Set the starting point to zero if it was not specified
  if (startingPoint == nullptr)
  {
    for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
    {
      this->bestPoint[iDOF] = this->currentPoint[iDOF] = (PrecisionType) 0.;
    }
  }
  else // otherwise copy the input
  {
    for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
    {
      this->bestPoint[iDOF] = this->currentPoint[iDOF] = startingPoint[iDOF];
    }
  }
  
  // Get the initial objective function value 
  // if it was not evaluated previously
  if (this->usePreEvaluatedEvaluatedObjectiveFunctionValue)
  {
    this->bestObjectiveFunctionValue =
      this->currentObjectiveFunctionValue =
      this->initialObjectiveFunctionValue =
      this->preEvaluatedObjectiveFunctionValue;
    // Reset to default behaviour
    this->usePreEvaluatedEvaluatedObjectiveFunctionValue = false;
  }
  else
  {
    this->bestObjectiveFunctionValue =
      this->currentObjectiveFunctionValue =
      this->initialObjectiveFunctionValue =
      this->objectiveFunction->GetValue( currentPoint );
  }

  // Print the objective function value at the beginning
  std::cout << "Optimiser --> it=[  0] fVal=" << std::setprecision( 5 ) << std::setw( 10 ) << this->bestObjectiveFunctionValue << std::endl;

  // Get the normalised gradient (i.e. fixed bool value)
  this->objectiveFunction->GetGradient( this->currentPoint, this->gradient, true );

  // Save the calculated gradient if required. 
  this->SaveIntermediateGradientImage();

  // Only set these parameters to the starting configuration if we do not 
  // want to continue a previous optimisation. Otherwise calcualte the 
  // conjugate gradient update. 
  if (continueOptimisation && this->continuationOfOptimisationPossible)
  {
    this->CalculateGradientUpdate();
  }
  else
  {
    for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
    {
      // Following reg-resp implementation here
      // numerical recipes would be  (for minimisation):
      //   array1[iDOF] = -gradient[iDOF]; 
      //   gradient[iDOF] = array2[iDOF] = array1[iDOF];
      this->array2[iDOF] = this->array1[iDOF] = -this->gradient[iDOF];
    }
    // Indicate that the arrays were initialised correctly and gradient updates can be used from here on
    this->continuationOfOptimisationPossible = true;
  }

  // define the step size
  PrecisionType currentLineSearchStepSize = this->maxLineSearchStepSize;

  for (this->currentIterationNumber = 1; (this->currentIterationNumber < this->maxNumberOfIterations) && (currentLineSearchStepSize > 0); this->currentIterationNumber++)
  {
    // Remember how many iterations have been performed
    // Note: Counter performs differently from the nifty-reg implementaiton
    this->numberOfIterationsPerformed = this->currentIterationNumber;

    // Perform the line maximisation along (Numerical Recipes: xi)
    // reg_resp: bestDOF[i] + scale * gradient[i]
    // note that nifty_reg_resp counted line optimisation iterations as iterations
    this->PerformLineSearch( this->maxLineSearchStepSize, this->minLineSearchStepSize, currentLineSearchStepSize );

    std::cout << "Optimiser --> it=[" << std::setw( 3 ) << this->currentIterationNumber << "]"
      << " fVal=" << std::setprecision( 5 ) << std::setw( 10 ) << this->bestObjectiveFunctionValue
      << " l=" << std::setprecision( 5 ) << std::setw( 10 ) << currentLineSearchStepSize << std::endl;

    // normal return check as definde in numerical recipes ignored in reg_resp
    // Note: Iteration number could have been reached here, hence check and 
    //       exit if necessary to prevent unneccessary calculation of gradient.
    if (this->currentIterationNumber >= this->maxNumberOfIterations) break;
    
    // No need to update the gradient if we don't use it thereafter
    if (currentLineSearchStepSize <= 0) break;

    // get the gradient at the current point (after line search this is eqaul to best point)
    this->objectiveFunction->GetGradient( this->currentPoint, this->gradient, true );
    
    // Save the calculated gradient if required. Normalise to make the output equivalent to the call before
    this->SaveIntermediateGradientImage();

    // Claculate the CG update
    this->CalculateGradientUpdate();
  }
}




//----------------------------------------
// ConjugateGradientOptimiser::Initialise
//----------------------------------------
void ConjugateGradientOptimiser::Initialise()
{
  if (this->objectiveFunction == nullptr)
  {
    supremo_print_error( "Objective function has to be set prior to initialisation of the optimiser." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }
  
  // Initialise step size of base class
  this->maxLineSearchStepSize = this->objectiveFunction->GetMaxStepSize();
  this->minLineSearchStepSize = this->maxLineSearchStepSize / (PrecisionType) 100.;

  this->numberOfDOFs = objectiveFunction->GetNumberOfParameters();
  
  this->array1       = (PrecisionType*)malloc( this->numberOfDOFs * sizeof( PrecisionType ) );
  this->array2       = (PrecisionType*)malloc( this->numberOfDOFs * sizeof( PrecisionType ) );
  this->gradient     = (PrecisionType*)malloc( this->numberOfDOFs * sizeof( PrecisionType ) );
  this->currentPoint = (PrecisionType*)malloc( this->numberOfDOFs * sizeof( PrecisionType ) );
  this->bestPoint    = (PrecisionType*)malloc( this->numberOfDOFs * sizeof( PrecisionType ) );
}




//-----------------------------------------------------
// ConjugateGradientOptimiser::CalculateGradientUpdate
//-----------------------------------------------------
void ConjugateGradientOptimiser::CalculateGradientUpdate()
{
  // Claculate the CG update as described in numerical recipes
  double gg = 0.0, dgg = 0.0;

  for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
  {
    gg += this->array2[iDOF] * this->array1[iDOF];
    dgg += (this->gradient[iDOF] + this->array1[iDOF]) * this->gradient[iDOF];
  }

  double gam = dgg / gg;

  for (unsigned int iDOF = 0; iDOF < this->numberOfDOFs; ++iDOF)
  {
    this->array1[iDOF] = -this->gradient[iDOF];
    this->array2[iDOF] = (this->array1[iDOF] + gam * this->array2[iDOF]);
    this->gradient[iDOF] = -this->array2[iDOF];
  }
}



void ConjugateGradientOptimiser::SaveIntermediateGradientImage()
{
  // Save the calculated gradient if required. Normalise to make the output equivalent to the call before
  if (!this->outputIntermediateGradientFolder.empty())
  {
    // Get the current gradient image(s) from the objective function
    std::vector<nifti_image*> intermediateGradImages = this->objectiveFunction->GetGradientAsImage( this->currentPoint, true );
    
    // Save each image individually
    for (size_t i = 0; i < intermediateGradImages.size(); ++i)
    {
      // The output file name. Note that we append an indentifier for each image according to the transformation
      std::ostringstream osOutFileName;
      osOutFileName << this->outputIntermediateGradientFolder
        << "/nmm_Grad_" << this->outputIntermediateGradientPrefix
        << "_t" <<std::setfill( '0' ) << std::setw( 2 ) << i 
        << "__" << std::setfill( '0' ) << std::setw( 2 )
        << this->currentIterationNumber << ".nii.gz";
      reg_io_WriteImageFile( intermediateGradImages[i], osOutFileName.str().c_str() );

      // We don't need the image after it was saved so it can be deleted. 
      nifti_image_free( intermediateGradImages[i] );
      intermediateGradImages[i] = nullptr;
    }
  }

  

}
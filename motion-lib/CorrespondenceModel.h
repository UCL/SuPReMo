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
#include <vector>
#include <iostream>
#include "Supremo.h"
#include "Transformation.h"


/** Class implementing the (respiratory) correspondence model.  
 */
class CorrespondenceModel{
public:

  typedef float PrecisionType;
  typedef std::vector<PrecisionType> SurrogateSignalType;
  
  /** Constructor.
   *  \param numberOfSurrogateSignalsIn The number of surrogate signals for this correspondence model.
   *  \param transformationIn The base transformation that is internally used to generate a transformation from input surrogate signals
   */
  CorrespondenceModel( unsigned int numberOfSurrogateSignalsIn, 
                       std::shared_ptr<Transformation> transformationIn );

  /** Destructor.
   */
  ~CorrespondenceModel();

  /** Calculate the transformation parameters for a given surrogate signal on the 
   *  basis of the current model parameters. 
   *  \f$ \mathbf{M}_t(\mathbf{S}_t) \f$.
   *  \param surrogateSignalIn Pointer to the surrogate values \f$ \mathbf{S}_t \f$.
   */
  //PrecisionType* GetTransformParamsFromSurrogateSignal( PrecisionType* surrogateSignalIn );

  /** Get a transformatinon object corresponding to a surrogate signal. Checks the size of the surrogate signal. 
   */
  std::shared_ptr<Transformation> GetTransformationFromSurrogateSignal( const SurrogateSignalType & surrogateSignalIn );

  /** Get access to the current parameters.
   *  A pointer to  \f$ \mathbf{R} \f$ is returned.
   */
  PrecisionType* GetParameters();

  /** Set the current parameters.
   *  The internal data of \f$ \mathbf{R} \f$ will be updated (i.e. copied).
   */
  void SetParameters( const PrecisionType* parametersIn );


  /** Set the current parameters from the given input images. 
   *  This function takes the image intensity and sets it as the correspondence model parameters. At the moment 
   *  only the total number of parameters is tested. 
   *  \param parameterImagesIn Vectore with images holding the correspondence model parameters.
   */
  void SetParameters( const std::vector<nifti_image*>& parameterImagesIn );


  /** Function used to calculate the gradient length, which in turn is used to normalise the gradient of the objective function. 
   *  The correspondence model class defines the relationship between the respiratory correspondence model parameters and the 
   *  transformation parameters. This is dependent on the number of surrogate signals. The maximum transformation parameter length
   *  will be calculated by the member transformation. 
   */
  PrecisionType GetMaxParameterLength( PrecisionType* parametersIn );

  /** Get the number of parameters.
   */
  inline unsigned int GetNumberOfParameters() { return this->numberOfModelParameters; };

  /** Get the gradient of the transformation parametres with respect to the model parametres.
   *  This function implements \f$ \frac{\partial \mathbf{M}_t }{\partial \mathbf{R} } \f$.
   *  Adds calculates the correspondence model gradient 
   */
  void GetTransformationParameterGradientWRTModelParameters( PrecisionType* transformationGradientIn,
                                                             const SurrogateSignalType & surrogateSignalIn,
                                                             PrecisionType* correspondenceModelGradientOut );

  /** Initialise the given level. 
   *  The transformation will be used to upsample the model parameters.
   */
  void InitialiseLevel( unsigned int levelIn );

  /** Get a pointer to the transformation. Allow other objects, such as the ObjectiveFunction to use
   *  methods provided by the transformation
   */
  std::shared_ptr<Transformation> GetTransformation();
  
  /** Generates a vector with nifti images which can be saved.
   *  The data of the images will be allocated for the images, so they are detached from 
   *  the internal representation of the correspondence model. 
   *  The transformations can have multiple images so each image returned by this function 
   *  reflects the number of images of the underlying transformation. 
   */
  std::vector<nifti_image*> GetCorrespondenceModelAsImage();

  /** Save the current model parameters to an internal variable such that they can be recovered if necessary
   *  by CorrespondenceModel::RecoverSavedModelParameters().
   */
  void SaveCurrentModelParametersForRecovery();


  /** Recover correspondence model parameters that were previously saved. Parameters for recovery have to be
   *  saved previously using CorrespondenceModel::SaveCurrentModelParametersForRecovery().
   */
  void RecoverSavedModelParameters();


private:
  /** Free the memory allocated to store the recovery model parameters
   *  and set the corresponding pointer to null.
   */
  void ClearRecoveryModelParameters();

  unsigned int numberOfSurrogateSignals;          ///< The number of surrogate signals used.
  unsigned int numberOfTransformationParameters;  ///< The number of parameters to describe the transformation completely
  unsigned int numberOfModelParameters;           ///< The total number of parameters describing the correspondence model
  std::shared_ptr<Transformation> transform;      ///< The internal transformation object that is used to construct a transformation (copy) from a surrogate signal
  PrecisionType* modelParameters;                 ///< Pointer to the model parameters
  PrecisionType* recoveryModelParameters;         ///< Pointer to model parameters that may be recovered
};

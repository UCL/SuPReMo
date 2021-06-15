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



//----------
// includes
//----------
#include "Supremo.h"
#include <memory>
#include "SupremoUtils.h"



//----------------------
// Forward declarations
//----------------------
class CorrespondenceModel;
class Transformation;




/** Class to compute the objective function. The weigths of the similarity and penatly terms are computed internally.
 */
class ObjectiveFunction{
public:

  typedef float PrecisionType;
  typedef std::vector< std::vector<PrecisionType> > SurrogateSignalType;

  /** Constructor.
   */
  ObjectiveFunction();


  /** Destructor.
   */
  ~ObjectiveFunction();


  /** Get the value of the objective function for a specific set of parameters
   *  Calculates \f$ \mathcal{C}_\text{total} = \sum_{t=1}^{N_i} \lambda \mathcal{S}( \mathbf{P}_t, A_t(\mathbf{T}(\mathbf{I}, \mathbf{M}_t)))
               + (1-\lambda)\mathcal{R}(\mathbf{M}_t) \f$
   *  \param parametersIn Pointer to the parameters that will be fed into the correspondence model, or point \f$\mathbf{R}\f$ in high-dimensional space for 
   *                      which the objective function-value will be calculated.
   */
  PrecisionType GetValue( const PrecisionType* parametersIn );


  /** Get the gradient of the objective function for a specific set of parameters
   *  This function copies the gradienbt of the objective function to an externally allocated memory location. 
   * 
   *  \param parametersIn Pointer to the parameters that will be fed into the correspondence model
   *  \param gradientOut  Pointer to where the gradient will be written. Calling function needs to handle allocation of gradient.
   *  \param normaliseGradient Set to true if the gradient should be normalised with the maximum transformation length
   */
  void GetGradient( const PrecisionType* parametersIn, PrecisionType* gradientOut, bool normaliseGradient = false );


  /** Calculate the gradient for the input parameters and return an image with the gradient
   *  \param parametersIn Pointer to the parameters that will be fed into the correspondence model
   *  \param normaliseGradient Set to true of the gradient should be normalised. Also see \ref GetGradient().
   */
  std::vector<nifti_image*> GetGradientAsImage( const PrecisionType* parametersIn, bool normaliseGradient = false );


  /** Get the number of parameters (degrees of freedom) of the objectve function.
   */
  unsigned int GetNumberOfParameters();


  /** Get the maximum step size for the line-optimisation. 
   *  Calcualtes the maximum step size \f$ l_\mathrm{max} \f$ from the voxel size of the reference-state image (\f$ \{d_x, d_y, d_z\} \f$) and the surrogate signal values \f$ s_i \f$. 
   *  \f[ l_\mathrm{max} = \frac{\max(\{d_x, d_y, d_z\})}{\max( \| s_i \| ) } \f]
   *  Obsiously both, surrogate signal and reference-state image have to be set before calling this function. 
   */
  PrecisionType GetMaxStepSize();


  /** Set the correspondence model. When the correspondence model is set, the transformation is defined and 
   *  the similarity and penalty term weights are updated.
   *  \param correspondenceModelIn Pointer to the correspondence model used
   */
  void SetCorrespondenceModel( const std::shared_ptr<CorrespondenceModel>& correspondenceModelIn );


  /** Set the object that measures the similarity between images.
   *  \param imageSimilarityIn Image similarity object. 
   */
  void SetSimilarityMeasure( const std::shared_ptr<ImageSimilarity>& imageSimilarityIn);


  /** Set the surrogate signal
   *  To measure the objective function, all the surrogate values need to be known (and of correct size).
   */
  void SetSurrogateSignals( const SurrogateSignalType & surrSignalsIn );


  /** Set the current reference state image
   *  \param refStateImgIn Pointer to the nifti_image structure holding the reference state image
   */
  void SetReferenceStateImage( nifti_image* refStateImgIn );


  /** Set the dynamic images. 
   *  Only set the resolution that is to be computed for each level. 
   *  \param dynamicImagesIn Vector holding the pointer to the dynamic images. Must be of same size as surrogate signals ObjectiveFunction::SetSurrogateSignals. 
   *  \param dynamicDataTypeIn Defines the type of the dynamic data
   */
  void SetDynamicImages( const std::vector<nifti_image*>& dynamicImagesIn, t_dynamicData dynamicDataTypeIn );

  /** Set the image acquisition object. Has to implement functionality to simulate the acquisition process and to calculate the 
   *  adjoint of that process. 
   *  \param imageAcquisitionIn Shared pointer to the image acquisition object.
   */
  void SetImageAcquisition( const std::shared_ptr<ImageAcquisition>& imageAcquisitionIn );

private:
  std::shared_ptr<CorrespondenceModel> correspondenceModel;   ///< Object that represents the correspondence model. Able to generate a transformation.
  std::shared_ptr<ImageSimilarity>  imageSimilarity;          ///< Object that measures the similarity between two images.
  SurrogateSignalType surrogateSignals;                       ///< All surrogate signals
  nifti_image* referenceStateImage;                           ///< Nifti image structure with all 
  std::vector<nifti_image*> dynamicImages;                    ///< Vector holding all pointers to the dynamic images
  t_dynamicData dynamicDataType;                              ///< Specifying dynamic data type
  PrecisionType similarityWeight;                             ///< Similarity weight \f[ 1-\lambda \f]
  std::shared_ptr<ImageAcquisition> imageAcquisition;         ///< The object simulating acquisition and calculating the adjoint of the image acquisition procedure.
};


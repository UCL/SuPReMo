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

#include "Supremo.h"
#include "nifti1_io.h"
#include "Transformation.h"





/** Class implementing the sliding transformation 
 */
class SlidingTransformation :public Transformation
{
public:
  typedef Transformation::PrecisionType PrecisionType;

  /** Constructor to generate a sliding transformation object. 
   *  
   *  The sliding transformation is defined by two individual transformation objects, i.e. transformationA and transformationB. 
   * 
   *  \param transformationAIn Pointer to the first transformation. 
   *  \param transformationBIn Pointer to the second transformation. 
   *  \param defSpaceImage Image defining the area where the sliding-transformation specific constraint term will be calculated. 
   *  \param numberOfLevels Number of pyramid levels. Required to correctly calculate the image pyramid of defSpaceImage.
   *  \param numberOfLevelsToPerform Number of levels to perform. See numberOfLevels.
   */
  SlidingTransformation( std::shared_ptr<Transformation>transformationAIn,
    std::shared_ptr<Transformation>transformationBIn,
    nifti_image* defSpaceImage, 
    unsigned int numberOfLevels, 
    unsigned int numberOfLevelsToPerform );
  
  /** Copy constructor.
   *  \param transformToCopy Reference to the sliding transformation object to copy.
   */
  SlidingTransformation( const SlidingTransformation& transformToCopy );

  /** Destructor.
   */
  ~SlidingTransformation();
  
  /** Initialise a given level. Call the initialisation for the individual transformations and update the current defSpaceImage. 
   *  \param level Defines the level that is to be initialised.
   */
  virtual void InitialiseLevel( unsigned int level );
  
  /** Set the parameters that define the sliding transformation.
   *  \param parametersIn Pointer to the parameters that will be copied to the internal transformation parameters.
   *  \param parametersAreDisplacements Set to true if the parameters are displacements. 
   */
  void SetParameters( PrecisionType* parametersIn, bool parametersAreDisplacements );
  
  /** For a given target image, calculate a deformation vector field (DVF). 
   *  \param targetImageIn The target image defines the geometry of the DVF, or the voxel locations for which deformation vectors are calculated.
   */
  virtual nifti_image* GetDeformationVectorField( nifti_image* targetImageIn );
   
  /** Calculates the gradient of an image for a given deformation vector field.
   *  \param denseDVFIn Deformation vector field 
   */
  virtual SlidingTransformation::PrecisionType* GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn );
  
  /** Calculate the gradient of the constraint term (regularisation) for the transformation with the current parameters.
   */
  virtual SlidingTransformation::PrecisionType* GetConstraintGradientWRTTransformationParameters();
  
  /** Get the constraint value for the transformation with the current parameters.
   */
  double GetConstraintValue();
  
  /** Return a pointer to the parameters of the transformation
   */
  virtual PrecisionType* GetCopyOfParameters();
  
  /** Function to create a deep copy of this object.
   */
  std::shared_ptr<Transformation> DeepCopy();
  
  /** Display function to print the transformation parameters to the standard output
   */
  void DisplayTransformationParameters();
  
  /** Calculates the sum of the penalty term weights. Responsibility of the correct 
   *  scaling of the penalty term lies with the generating function.
   */
  virtual PrecisionType GetSumOfPenaltyWeights();

  /** Calculate the maximum control-point displacement value for the given transformation parametes. 
   */
  virtual PrecisionType GetMaxTransformationParameterLength( PrecisionType* parametersIn );


  /** Returns the transformation as an image. 
   */
  std::vector<nifti_image*> GetTransformationAsImage();

  /** Set the sliding-transformation specific constraint weight. Higher values will result in less gaps and overlaps. 
   * 
   *  \param gapOverlabConstraintWeightIn Weight of the gap-overlap constraint term.
   */
  void SetGapOverlapConstraintWeight( double gapOverlabConstraintWeightIn );

  /** Set the signed distance map image used to determine the sliding boundary. This function generates a local copy of this image
   *  \param signedDistanceMapImageIn  Nifti image of the signed distance map in the source image space.
   */
  void SetSignedDistanceMapImage( nifti_image* signedDistanceMapImageIn);


private:
  
  /** Delete this objet's copy of the signed distance image used.
   */
  void ClearSignedDistanceMap();

  /** Delete the transfomred signed distance maps
   */
  void ClearTransformedSignedDistanceMaps();
  
  /** Allocate the transfomred signed distance maps
   *  \param targetImage The target image defines the size of the allocated signed distance maps. 
   */
  void AllocateTransformedSignedDistanceMaps( nifti_image* targetImage );


  /** Calculate the gap and overlap constraint term from the transfomred signed distance maps. 
   *  \note The function GetDeformationVectorField() has to be called before this can be computed.  
   */
  double CalculateGapOverlapConstraintTerm();

  /** Checks the geometry of the transformed signed distance map and the (current) image defining the 
   *  space where the sliding constraint is calculated. 
   */
  bool TransformedSignedDistMapsCanBeUsedForConstraintCalculations();
  
  std::shared_ptr<Transformation> transformA;   ///< The transformation defining the inside, i.e. SDT < 0
  std::shared_ptr<Transformation> transformB;   ///< The transformation defining the outside, i.e. SDT > 0
  
  nifti_image* signedDistMapImage;              ///< Nifti image holding the singed distance image in source-image space for the current level.
  
  std::shared_ptr<ImagePyramid<PrecisionType>> defSpaceImagePyramid; ///< The image pyramid holding the images in which the transformation is defined. Should correspond with those used to initilase the transfomrations A and B. Will be generated during construction. 
  nifti_image* curDefSpaceImage;                ///< Pointer to the current image within the defSpaceImagePyramid.  

  int lastInitialisedLevel;                     ///< The last initialised level. This determines which level of the def space image pyramid is accessed. Will be set to -1 before initialisation.

  nifti_image* signedDistMapImageTransformedA;  ///< Nfiti image holdeing the signed distance map transformed with transformA
  nifti_image* signedDistMapImageTransformedB;  ///< Nfiti image holdeing the signed distance map transformed with transformA
  
  double gapOverlapConstraintWeight;            ///< The gap-overlap-constraint weight
};


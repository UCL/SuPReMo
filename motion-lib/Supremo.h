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
#include "nifti1_io.h"
#include "ImagePyramid.h"
#include "ImageSimilarity.h"
#include "SSDImageSimilarity.h"
#include "CorrespondenceModel.h"
#include "Transformation.h"
#include "BSplineTransformation.h"
#include "SlidingTransformation.h"
#include "MoCoRecon.h"
#include "ObjectiveFunction.h"
#include "Optimiser.h"
#include "ConjugateGradientOptimiser.h"
#include "ImageAcquisition.h"
#include "NoImageAcquisition.h"
#include "LowResolutionImageAcquisition.h"


// Forward declaration
class CorrespondenceModel;

inline void supremo_exit(int Val, const char* fileName, int lineNumber)
{
  fprintf(stderr, "[Supremo] Exit here. File: %s:%i\n", fileName, lineNumber);
  exit(Val);
};

inline void supremo_print_debug(const char* message){ printf("[Supremo DEBUG] Message: %s\n", message); };
inline void supremo_print_warning(const char* message){ printf("[Supremo WARNING] Message: %s\n", message); };
inline void supremo_print_error(const char* message){ fprintf(stderr, "[Supremo ERROR] Message: %s\n", message); };



/** Class that implements the generalised motion modelling framework by McClelland et al. (2017).
 *  This class implements the generalised motion modelling framework and is the main class that determines flow control, data handling etc. 
 */
class Supremo
{
public:
  typedef float VoxelType;
  
  /** Constructor */
  Supremo();

  /** Destructor */
  ~Supremo();

  /** Fit a motion model and perform motion compensated image reconstruction */ 
  void FitMotionModelAndReconstruct();


  /** Set all dynamic images. 
   *  Used to set a number of externally loaded dynamic images.
   *  \param dynamicImagesIn Pointer to the dynamic images. 
   *  \param numberOfDynamicImagesIn Number of dynamic images, i.e. length of the pointer above.  
   */
  void SetDynamicImages(nifti_image** dynamicImagesIn, int numberOfDynamicImagesIn);


  /** Set the reference state image \f$ \mathbf{I}_0 \f$. This image was known as static or source image. 
   *  \param referenceStateImageIn Pointer to the nifti_image holding the static image. 
   */
  void SetReferenceStateImage(nifti_image* referenceStateImageIn);


  /** Set the surrogate data. 
   *  \note This function must be called after the dynamic images were set. 
   *  \param surrogateSignalsIn Pointer to the surrogate data. Has to be of size numberOfSurrogateSignalsIn*numberOfDynamicImages.
   *  \param numberOfSurrogateSignalsIn Number of surrogate signals per dynamic image.
   */
  void SetSurrogateSignals( float* surrogateSignalsIn, int numberOfSurrogateSignalsIn );


  /** Set an image to define the space of the deformed images. 
   *  Image used to define the space of the deformed images. This image is used to define the extent and resolution 
   *  (if specified in voxels) of the model and transform CPGs. If a defSpace image is not specified the reference state image 
   *  will be used.
   *  \param defSpaceImageIn Pointer to a nifti image that defines the space of the deforemd images.
   */
  void SetDefSpaceImage( nifti_image* defSpaceImageIn );


  /** Set the dynamic image type, i.e. same or lower resolution as static image
   *  \param dynamicImgaeDataTypeIn Type of the dynamic data.
   */
  void SetDynamicImageDataType( t_dynamicData dynamicImgaeDataTypeIn );

  /** Set the motion compensated image reconstruction class.
   *  \param mocoReconstructionIn Object that implements a motion compensated image reconstruction.
   */
  void SetMotionCompensatedReconstruction( std::shared_ptr<MoCoRecon> mocoReconstructionIn );

  /** Set the folder to which intermediate MCR results will be saved
   *  \param outputFolder
   */
  void SetInterMCROutputFolder( const std::string& outputFolder );

  /** Set the folder to which intermediate objective function gradients will be saved
     *  \param outputFolder
     */
  void SetInterGradOutputFolder( const std::string& outputFolder );

  /** Set the transformation type used
   *  \param transformationTypeIn Transformation type
   */
  void SetTransformationType( t_transformationType transformationTypeIn );

  /** Set the distance map for the sliding transformation
   *  \param slidingTrafoDistanceMapIn Pointer to the nifti image structure holding the sliding transformation
   */
  void SetSlidingTrafoDistanceMap( nifti_image* slidingTrafoDistanceMapIn );


  /** Set the gap/overlap constraint weight. 
   *  This weight only has an effect if the sliding transformation is being used. 
   *  \param gapOverlapConstraintWeightIn Gap/overlap constraint weight. 
   */
  void SetSlidingGapOverlapConstraintWeight( float gapOverlapConstraintWeightIn );


  /** Set the control point spacing  of the B-spline control point grid
   *  \param sx Control point spacing in x direction
   *  \param sy Control point spacing in y direction
   *  \param sz Control point spacing in z direction
   */
  void SetBSplineCPGSpacing( float sx, float sy, float sz );

  /** Set the bending energy weight for the B-spline regularisation
   *  \param be Bending energy weight to set
   */
  void SetBSplineBendingEnergy( float be );

  /** Set the linear energy weight for the B-spline regularisation
   *  \param le Linear energy weight to set
   */
  void SetBSplineLinearEnergy( float le );

  /** Set the number of pyramid levels
   *  \param numberOfLevelsIn Total number of pyramid levels
   */
  inline void SetNumberOfPyramidLevels( unsigned int numberOfLevelsIn )
  {
    this->numberOfLevels = numberOfLevelsIn;
  }

  /** Set the number of pyramid levels on which a computation should be performed.
  *  \param numberOfLevelsToPerformIn Nmber of pyramid levels. Has to be equal or smaller than the value set with \ref SetNumberOfPyramidLevels(). 
  */
  inline void SetNumberOfPyramidLevelsToPerform( unsigned int numberOfLevelsToPerformIn )
  {
    this->numberOfLevelsToPerform = numberOfLevelsToPerformIn;
  }

  /** Set the number of iterations used to fit the motion model.
   *  \param maxModelFittingIterations Maximum number of model fitting iterations
   */
  inline void SetMaxModelFittingIterationNumber( unsigned int maxModelFittingIterations )
  {
    this->maxModelFitIterationNumber = maxModelFittingIterations;
  }


  /** Set the number of switch iterations. 
   *  Determines how often the method iterates between model fitting and motion-compensated image reconstruction. 
   *  \param maxSwitchIterations Maximum number of switch iterations
   */
  inline void SetMaxSwitchIterationNumber( unsigned int maxSwitchIterations )
  {
    this->maxSwitchIterationNumber = maxSwitchIterations;
  }

  /** Set the input RCM correspondence model. 
   *  The values of this image/these image(s) will be used as a starting point of the optimisation. Multiple RCMs can occur if the sliding 
   *  transforation was used. 
   *  \param inputRCMImgs Vector holding image pointers to the input respiratory correspondence model images.
   */
  inline void SetInputRCMImages(const std::vector<nifti_image*>& inputRCMImgs)
  {
    this->inputRCMImages = inputRCMImgs;
  };



  /** Get the respiratory correspondence model after model fitting 
   *   \return Vecotr with pointers to the nifti images representing the correspondence model. These will usually be a 5d nifti images.
   *   \note   The correspondence model image depends on the transformation model used. 
   */
  std::vector<nifti_image*> GetCorrespondenceModelAsImage();


  /** Get the resulting motion-compensated reconstructed image. 
   *  \return Pointer to the image holding the final MCR iamge. This image will be deallocated by this class on destruction. 
   */
  nifti_image* GetMotionCompensatedReconstructedImage() { return this->finalMoCoReconImage; }

  
  /** Once a motion model was fitted (or input), use the results to simulate the dynamic images
   *  \return A vector holding the pointers to the simulated dynamic images
   */
  std::vector<nifti_image*> SimulateDynamicImages();


  /** Once a motion model was fitted (or input), use the results to calculate the deformation vector fields. 
   *  The requested DVFs will have the size of the dynamic images.
   *  \return A vector holding the pointers to the DVFs for the defined surrogate values.
   */
  std::vector<nifti_image*> GenerateDVFsFromCorrespondenceModel();


private:
  /** Initillise the fitting and reconstruction method
   *   - Checks that all parameters were set correctly. 
   *   - Generates the image pyramids.
   */
  void Initialise(); 


  /** Checks that all parameters for fitting and reconstruction were provided correctly.
   */
  void CheckParameters();

  /** Display some parameters of the current levels including image resolution, transformation parameters, etc. 
   *  \param level The level number
   */
  void DisplayCurrentLevelParameters( unsigned int level );


  /** Save the current MCR image to file. The string levelID info will be used to determine the output
   *  \param mcrImage The reconstructed image that will be saved to the defined output directory. 
   */
  void SaveInterMCRToImage(nifti_image* mcrImage);

  unsigned int numberOfSurrogateSignals;          ///< Number of surrogate signals (per dynamic image)
  unsigned int numberOfDynamicImages;             ///< Total number of dynamic images
  
  nifti_image ** allDynamicImages;                ///< Pointer to all dynamic nifti images available
  nifti_image * referenceStateImage;              ///< Pointer to the reference state image \f$ \mathbf{I}_0 \f$
  nifti_image * finalMoCoReconImage;              ///< Pointer holding the final motion-compensated reconstruction image after fit and reconstruct 
  
  std::vector<std::vector<float> > surrogateSignals; ///< Vector of vector containing the the surrogate signals. First index for time point, second for surrogate signal number \verbatim surrogateSignals[time][surrNumber] \endverbatim
  
  nifti_image* defSpaceImage;                     ///< Pointer to image that defines the space of the deformed images.

  t_dynamicData dynamicImageDataType;             ///< Dynamic data type 

  std::vector<nifti_image*> inputRCMImages;       ///< The vector holding the pointer to the respiratory correspondence model images. Will be used as a starting point of the optimisation. 

  unsigned int maxSwitchIterationNumber;          ///< Maximum number of times to iterate between motion compensated reconstruction and fitting the respiratory correspondence model
  unsigned int maxModelFitIterationNumber;        ///< Maximum number of respiratory correspondence model fitting iterations

  unsigned int numberOfLevels;                    ///< The total number of levels
  unsigned int numberOfLevelsToPerform;           ///< The number of levels that are used in the fitting/reconstruction process
  unsigned int currentLevel;                      ///< The current level that is being processed

  t_transformationType transformationType;        ///< The transformation type used to deform the reference-state image

  float bSplineCPGSpacing[3];                     ///< B-spline control point spacing (used if b-spline transformation is used)
  float bSplineBendingEnergyWeight;               ///< Bending energy weight to constrain a b-spline transformation
  float bSplineLinearEnergyWeight;                ///< Linear energy weight to constrain a b-spline transformation

  nifti_image* slidingSignedDistanceMapImage;     ///< Image required by the sliding transformation to determine the regions involved
  float slidingGapOverlapConstraintWeight;        ///< The weight with which the gaps and overlaps will be constrained

  std::string outputInterMCRFolder;               ///< Folder to which intermediate MCRs will be saved
  std::string outputInterGradientFolder;          ///< Folder to which intermediate objective function gradients will be saved

  bool initialised;                               ///< Control parameter indicating if the fit-and-reconstruct method was initialised.

  // Image pyramids
  std::shared_ptr<ImagePyramid<VoxelType> > referenceStatePyramid;            ///< Pointer holding the reference state image pyramid
  std::vector<std::shared_ptr<ImagePyramid<VoxelType> > > allDynamicPyramids; ///< Vector of pointers holding all the dynamic image pyramids

  // Temporary image pointers
  nifti_image* currentReferenceStateImage;

  std::shared_ptr<Transformation> transform;                  ///< Pointer to the transformation used (e.g. b-spline)
  std::shared_ptr<CorrespondenceModel> correspondenceModel;   ///< Pointer to the correspondence model
  std::shared_ptr<ImageSimilarity> similarityMeasure;         ///< Pointer to the image similarity measure used
  std::shared_ptr<ImageAcquisition> imageAcquisition;         ///< Pointer to the image acquisition used
  std::shared_ptr<MoCoRecon> mocoReconstructor;               ///< The motion-compensated image reconstruction object

  std::string levelID;                                        ///< String holding an ID to differentiate file output from various multi-resolution levels and MCR iterations
};

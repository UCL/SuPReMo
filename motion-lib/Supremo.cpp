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
#include "Supremo.h"
#include "_reg_ReadWriteImage.h"
#include "_reg_tools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>  




//------------------
// Supremo::Supremo
//------------------
Supremo::Supremo() :
  numberOfSurrogateSignals( 0 ),
  numberOfDynamicImages( 0 ),
  allDynamicImages( nullptr ),
  referenceStateImage( nullptr ),
  defSpaceImage( nullptr ),
  finalMoCoReconImage( nullptr ),
  dynamicImageDataType( SAME_RES_AS_STATIC ),
  maxSwitchIterationNumber( -1 ),
  maxModelFitIterationNumber( 0 ),
  numberOfLevels( 3 ),
  numberOfLevelsToPerform( 3 ),
  currentLevel( 0 ),
  transformationType( STANDARD_B_SPLINE ),
  bSplineBendingEnergyWeight( 0 ),
  bSplineLinearEnergyWeight( 0 ),
  slidingSignedDistanceMapImage( nullptr ),
  slidingGapOverlapConstraintWeight( 0 ),
  outputInterMCRFolder( "" ),
  outputInterGradientFolder( "" ),
  initialised( false ),
  referenceStatePyramid( nullptr ),
  currentReferenceStateImage( nullptr ),
  transform( nullptr ),
  correspondenceModel( nullptr ),
  similarityMeasure( nullptr ),
  imageAcquisition( nullptr ),
  mocoReconstructor( nullptr ),
  levelID( "" )
{
  // Array needs to be initialised here for compatibility with VS2013
  // otherwise bSplineSpacing{0,0,0} in list above would have worked.
  this->bSplineCPGSpacing[0] = this->bSplineCPGSpacing[1] = this->bSplineCPGSpacing[2] = 0;
}




//-------------------
// Supremo::~Supremo
//-------------------
Supremo::~Supremo()
{
  this->transform.reset();
  nifti_image_free( this->finalMoCoReconImage );
}




//---------------------------
// Supremo::SetDynamicImages
//---------------------------
void Supremo::SetDynamicImages(nifti_image** dynamicImagesIn, int numberOfDynamicImagesIn)
{
  // Save the dynamic images to the member variables
  this->allDynamicImages = dynamicImagesIn;
  this->numberOfDynamicImages = numberOfDynamicImagesIn;
  }




//---------------------------------
// Supremo::SetReferenceStateImage
//---------------------------------
void Supremo::SetReferenceStateImage( nifti_image* referenceStateImageIn )
{
  // Save the pointer to the reference state image to the member variable
  this->referenceStateImage = referenceStateImageIn;
  return;
}



//-----------------------------------
// Supremo::SetReferenceStateImage
//-----------------------------------
void Supremo::SetDefSpaceImage( nifti_image* defSpaceImageIn )
{
  // Save an image that defines the reference space
  this->defSpaceImage = defSpaceImageIn;
}



//--------------------------------
// Supremo::SetSurrogateSignals
//--------------------------------
void Supremo::SetSurrogateSignals( float* surrogateSignalsIn, int numberOfSurrogateSignalsIn )
{
  // Remember the number of surrogate signals
  this->numberOfSurrogateSignals = numberOfSurrogateSignalsIn;

  if (this->numberOfDynamicImages == 0)
  {
    char msg[200];
    sprintf_s( msg, "Number of dynamic images is %i. These have to be set prior to calling this function", this->numberOfDynamicImages );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Save the raw surrogate data pointer into an array
  for (int iTimePoint = 0; iTimePoint < this->numberOfDynamicImages; ++iTimePoint)
  {
    std::vector<float> tmpVect;
    for (int iSurrSig = 0; iSurrSig < this->numberOfSurrogateSignals; ++iSurrSig)
    {
      tmpVect.push_back( surrogateSignalsIn[iTimePoint*this->numberOfSurrogateSignals + iSurrSig] );
    }
    this->surrogateSignals.push_back( tmpVect );
  }
}



//------------------------------------
// Supremo::SetDynamicImageDataType
//------------------------------------
void Supremo::SetDynamicImageDataType( t_dynamicData dynamicImgaeDataTypeIn )
{
  // Set the dynamic image data type
  this->dynamicImageDataType = dynamicImgaeDataTypeIn;
  return;
}



//-----------------------
// Supremo::SetMCRType
//-----------------------
void Supremo::SetMotionCompensatedReconstruction( std::shared_ptr<MoCoRecon> mocoReconstructionIn )
{
  // Set the motion-compensated image reconstruction type
  this->mocoReconstructor = mocoReconstructionIn;
}




//------------------------------------
// Supremo::SetInterMCROutputFolder
//------------------------------------
void Supremo::SetInterMCROutputFolder( const std::string & outputFolderIn)
{
  // Define the folder name where to save the intermediate motion-compensated images
  this->outputInterMCRFolder = outputFolderIn;
}



//------------------------------------
// Supremo::SetInterMCROutputFolder
//------------------------------------
void Supremo::SetInterGradOutputFolder( const std::string & outputFolderIn )
{
  // Define the folder name where to save the intermediate gradient images will be saved
  this->outputInterGradientFolder = outputFolderIn;
}




//------------------------------------
// Supremo::SetInterMCROutputFolder
//------------------------------------
void Supremo::SetTransformationType( t_transformationType tranformationTypeIn )
{
  // Set the transformation type
  this->transformationType = tranformationTypeIn;
}



//---------------------------------
// Supremo::SetBSplineCPGSpacing
//---------------------------------
void Supremo::SetBSplineCPGSpacing( float sxIn, float syIn, float szIn )
{
  // Define the bspline spacing in x, y, and z
  this->bSplineCPGSpacing[0] = sxIn;
  this->bSplineCPGSpacing[1] = syIn;
  this->bSplineCPGSpacing[2] = szIn;
}



//----------------------------------
// Supremo::SetBSplineBendingEnergy
//----------------------------------
void Supremo::SetBSplineBendingEnergy( float beIn )
{
  // Define the bending energy (on the final level)
  this->bSplineBendingEnergyWeight = beIn;
}




//---------------------------------
// Supremo::SetBSplineLinearEnergy
//---------------------------------
void Supremo::SetBSplineLinearEnergy( float leIn )
{
  // Define the linear energy for the bspline transformation regularisation
  this->bSplineLinearEnergyWeight = leIn;
}




//-----------------------------------------------
// Supremo::SetSlidingGapOverlapConstraintWeight
//-----------------------------------------------
void Supremo::SetSlidingGapOverlapConstraintWeight( float goIn )
{
  // Set the gap/overlap constraint weight for the sliding transformation
  this->slidingGapOverlapConstraintWeight = goIn;
}




//-------------------------------------
// Supremo::SetSlidingTrafoDistanceMap
//-------------------------------------
void Supremo::SetSlidingTrafoDistanceMap( nifti_image* signedDistanceMapImageIn )
{
  this->slidingSignedDistanceMapImage = signedDistanceMapImageIn;
}




//-----------------------------------------
// Supremo::FitMotionModelAndReconstruct
//-----------------------------------------
void Supremo::FitMotionModelAndReconstruct()
{
  // Initialise
  this->Initialise();
  
  // Loop over resolution levels
  for (; currentLevel < this->numberOfLevelsToPerform; ++currentLevel)
  {
    // Update the level ID to have the same ID for various outputs
    {
      std::ostringstream ossLevelID;
      ossLevelID.clear();
      ossLevelID << "lev_" << std::setfill( '0' ) << std::setw( 2 ) << currentLevel << "_initial";
      this->levelID = ossLevelID.str();
    }

    // Initialise the current level of the transformation and correspondence model
    this->transform->InitialiseLevel( currentLevel );
    this->correspondenceModel->InitialiseLevel( currentLevel );

    // member variable currentReferenceStateImage may not be required
    this->currentReferenceStateImage = this->referenceStatePyramid->GetLevel( currentLevel );

    // Display current level parameters
    this->DisplayCurrentLevelParameters( currentLevel );

    // Extract relevant dynamic images from  vector of pyramids
    std::vector<nifti_image*> curDynamicImages;
    for (int nDynImg = 0; nDynImg < this->numberOfDynamicImages; ++nDynImg)
    {
      curDynamicImages.push_back( this->allDynamicPyramids[nDynImg]->GetLevel( currentLevel ) );
    }

    // Motion compensated image reconstruction 
    if (this->mocoReconstructor != nullptr)
    {
      this->mocoReconstructor->SetSurrogateSignals( this->surrogateSignals );
      this->mocoReconstructor->SetDynamicImages( curDynamicImages );
      this->mocoReconstructor->SetCorrespondenceModel( this->correspondenceModel );
      this->mocoReconstructor->SetImageAcquisition( this->imageAcquisition );
      this->mocoReconstructor->SetReconstructionGeometryImage( this->currentReferenceStateImage );
      
      // and perform the first reconstruction for this level
      this->mocoReconstructor->Update();

      // copy the reconstructed image into the current reference state image
      this->mocoReconstructor->CopyReconstructedImageContentsToImage( this->currentReferenceStateImage );
      
      // Update the levelID used to save the intermediate results
      this->SaveInterMCRToImage( this->mocoReconstructor->GetReconstructedImage() );
    }

    // Initialise the objective function to be optimised and 
    // feed in all required paramters
    std::shared_ptr<ObjectiveFunction> objectiveFunction = std::make_shared<ObjectiveFunction>();
    objectiveFunction->SetCorrespondenceModel( this->correspondenceModel );
    objectiveFunction->SetSurrogateSignals( this->surrogateSignals );
    objectiveFunction->SetSimilarityMeasure( this->similarityMeasure );
    objectiveFunction->SetReferenceStateImage( this->currentReferenceStateImage );
    objectiveFunction->SetDynamicImages( curDynamicImages, this->dynamicImageDataType );
    objectiveFunction->SetImageAcquisition( this->imageAcquisition );

    // initialise the optimiser for fitting the model
    std::shared_ptr<Optimiser> optimiser = std::make_shared<ConjugateGradientOptimiser>();
    optimiser->SetObjectiveFunction( objectiveFunction );
    optimiser->SetMaxIterations( this->maxModelFitIterationNumber );
    optimiser->SetOutputIntermediateGradientFolder( outputInterGradientFolder, levelID );
    optimiser->Initialise();

    // Need to store objective function value from previous switch iteration (i.e. starting point of the optimisation)
    // Hence this will be read from the optimisation afterwards;
    double prevObjectiveFunctionValue;

    // Ask the correspondence model to save the the model parameters 
    // in case they need to be recovered later on
    this->correspondenceModel->SaveCurrentModelParametersForRecovery();
    
    // switch between model fitting and motion compensated image reconstruction
    // until max number of switch iterations performed
    for (int switchIter = 0; switchIter < this->maxSwitchIterationNumber; switchIter++)
    {
      // Update the level ID to have the same ID for various outputs
      {
        std::ostringstream ossLevelID;
        ossLevelID.clear();
        ossLevelID << "lev_" << std::setfill( '0' ) << std::setw( 2 ) << currentLevel
          << "_iter_" << std::setfill( '0' ) << std::setw( 2 ) << switchIter;
        this->levelID = ossLevelID.str();
      }

      optimiser->SetOutputIntermediateGradientFolder( outputInterGradientFolder, levelID );
      
      // Perfrom the optimisation with the starting point defined by the correspondence model, 
      // iteration counts are reset everytime Optimise() is being called
      optimiser->Optimise( this->correspondenceModel->GetParameters() );
      this->correspondenceModel->SetParameters( optimiser->GetBestPoint() );

      prevObjectiveFunctionValue = optimiser->GetInitialObjectiveFunctionValue();

      //check if objective function improved from model fitting
      if (prevObjectiveFunctionValue >= optimiser->GetBestObjectiveFunctionValue())
      {
        // model has not changed (or is worse than before, but this should be impossible)
        // so no need to redo reconstruction unless super-resolution with update is requried
        if (nullptr == this->mocoReconstructor || !this->mocoReconstructor->GetRepeatedUpdateChangesMCR())
          break;
      }

      CorrespondenceModel::PrecisionType objectiveFunctionValueAfterRecon = optimiser->GetBestObjectiveFunctionValue();

      if (nullptr != mocoReconstructor)
      {
        // Update the motion-compensated image reconstruction
        // The correspondence model was already updated
        mocoReconstructor->Update();

        // Save the newly reconstructed image if requried
        this->SaveInterMCRToImage( mocoReconstructor->GetReconstructedImage() );

        // the updated reference state image will have an impact on the objective function value
        // ideally this has improved
        // need to evaluate the objective function 
        // but this is done when the optimiser is called... so if we do this here explicitly, 
        objectiveFunction->SetReferenceStateImage( mocoReconstructor->GetReconstructedImage() );
        objectiveFunctionValueAfterRecon = objectiveFunction->GetValue( this->correspondenceModel->GetParameters() );
        objectiveFunction->SetReferenceStateImage( this->currentReferenceStateImage );

        std::cout << "fValAfterReco=" << std::setprecision( 5 ) << objectiveFunctionValueAfterRecon << std::endl;
      }

      std::cout << "End of switch iteration " << switchIter << std::endl;
      
      // if the reconstruction results in a better objective function than the one directly after model fitting, then 
      // update the correspondence model and the reference state image
      if (objectiveFunctionValueAfterRecon > prevObjectiveFunctionValue) 
      {
        // updating the reconstruction improved the objective function, so update and continue
        if (nullptr != mocoReconstructor) mocoReconstructor->CopyReconstructedImageContentsToImage( this->currentReferenceStateImage );
        
        this->correspondenceModel->SaveCurrentModelParametersForRecovery();
        prevObjectiveFunctionValue = objectiveFunctionValueAfterRecon;
        optimiser->SetPreviouslyEvaluatedObjectiveFunctionValue( objectiveFunctionValueAfterRecon );
      }
      else
      {
        // The fitting and reconstruction did not improve the result. Do NOT upadte the current 
        // reference-state image and recover the correspondence model parameters
        this->correspondenceModel->RecoverSavedModelParameters();

        break;
      }
    } // for (int switchIter = 0; switchIter < this->maxSwitchIterationNumber; switchIter++)
  } // for currentLevel

  // The current reference state image is located in the image pyramid. Hence, we better copy it into a dedicated 
  // destination instead. 
  this->finalMoCoReconImage = nifti_copy_nim_info( this->currentReferenceStateImage );
  this->finalMoCoReconImage->data = malloc( this->finalMoCoReconImage->nvox * this->finalMoCoReconImage->nbyper );
  memcpy( this->finalMoCoReconImage->data, 
    this->currentReferenceStateImage->data, 
    this->finalMoCoReconImage->nvox * this->finalMoCoReconImage->nbyper );

  return;
}




//----------------------------------------
// Supremo::GetCorrespondenceModelAsImage
//----------------------------------------
std::vector<nifti_image*> Supremo::GetCorrespondenceModelAsImage()
{
  // Cannot get the correspondence model if it does not exist (anymore)
  if (this->correspondenceModel == nullptr)
  {
    supremo_print_error( "No correspondence model set." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Request the image from the correspondence model
  std::vector<nifti_image*> correspModelImages = this->correspondenceModel->GetCorrespondenceModelAsImage();
  
  // Set description in the nifti header accordingly in each returned image
  for (size_t i = 0; i < correspModelImages.size(); ++i)
  {
    memset( correspModelImages[i]->descrip, 0, 80 );
    strcpy( correspModelImages[i]->descrip, "Respiratory correspondence model from SuPReMo" );
  }

  return correspModelImages;
}




//--------------------------------
// Supremo::SimulateDynamicImages
//--------------------------------
std::vector<nifti_image*> Supremo::SimulateDynamicImages()
{
  
  // Check that all required members are available
  if (nullptr == this->allDynamicImages)
  {
    supremo_print_error( "Dynamic images not available." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  
  if (nullptr == this->imageAcquisition)
  {
    supremo_print_error( "Image acquisition not defined." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  if (nullptr == this->correspondenceModel)
  {
    supremo_print_error( "Correspondence model not defined." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  if (nullptr == this->finalMoCoReconImage)
  {
    supremo_print_error( "Final motion-compensated reconstructed image not defined." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  
  std::vector<nifti_image*> outDynamicImages;

  for (unsigned int iImg = 0; iImg < this->numberOfDynamicImages; ++iImg)
  {
    // Get the transformation for the current surrogate value
    std::shared_ptr<Transformation> curTrafo = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[iImg] );

    // Convert the current dynamic image to the data type needed
    nifti_image* curDynImg = nifti_copy_nim_info( this->allDynamicImages[iImg] );
    curDynImg->data = malloc( curDynImg->nbyper * curDynImg->nvox );
    memcpy( curDynImg->data, this->allDynamicImages[iImg]->data, curDynImg->nbyper * curDynImg->nvox );
    reg_tools_changeDatatype<VoxelType>( curDynImg );

    // Use the image acquisition class to define the minimum sized image that we need to define the 
    nifti_image* curMinSizeImgInFullImgSpace = this->imageAcquisition->AllocateMinimumSizeImgInFullImgSpace( this->finalMoCoReconImage, curDynImg, iImg );
    
    // Use the transformation to warp the final image into the min-sized space
    curTrafo->TransformImage( this->finalMoCoReconImage, curMinSizeImgInFullImgSpace );

    // Simulate the acquisition
    nifti_image* curSimulatedDynamicImg = this->imageAcquisition->SimulateImageAcquisition( curMinSizeImgInFullImgSpace, curDynImg, iImg);
    outDynamicImages.push_back( curSimulatedDynamicImg );

    // Free the intermediate images if it was not simply copied over by the image-acquisition object
    if (curMinSizeImgInFullImgSpace != curSimulatedDynamicImg)
    {
      nifti_image_free( curMinSizeImgInFullImgSpace );
    }
    nifti_image_free( curDynImg );
  }

  return outDynamicImages;
}



//----------------------------------------------
// Supremo::GenerateDVFsFromCorrespondenceModel
//----------------------------------------------
std::vector<nifti_image*> Supremo::GenerateDVFsFromCorrespondenceModel()
{

  // Check that all required members are available
  if (nullptr == this->allDynamicImages)
  {
    supremo_print_error( "Dynamic images not available." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  if (nullptr == this->imageAcquisition)
  {
    supremo_print_error( "Image acquisition not defined." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  if (nullptr == this->correspondenceModel)
  {
    supremo_print_error( "Correspondence model not defined." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Generate a vector that holds the generated DVF images
  std::vector<nifti_image*> outDVFs;

  for (unsigned int iImg = 0; iImg < this->numberOfDynamicImages; ++iImg)
  {
    // Get the transformation for the current surrogate value
    std::shared_ptr<Transformation> curTrafo = this->correspondenceModel->GetTransformationFromSurrogateSignal( this->surrogateSignals[iImg] );

    // The transformation will return a pointer to the internal nifti image holding the dvf. 
    // Request it for the size of the current dynamic image
    // Hence  a copy of this DVF needs to be generated, otherwise it will be lost
    nifti_image* curDVF = curTrafo->GetDeformationVectorField( this->allDynamicImages[iImg] );

    nifti_image* curDVFOut = nifti_copy_nim_info( curDVF );
    curDVFOut->data = malloc( curDVFOut->nvox * curDVFOut->nbyper );
    memcpy( curDVFOut->data, curDVF->data, curDVFOut->nvox * curDVFOut->nbyper );
    outDVFs.push_back( curDVFOut );
  }

  return outDVFs;
}




//---------------------
// Supremo::Initialise
//---------------------
void Supremo::Initialise()
{
  // No need to initialise if this was already done before
  if (this->initialised) return;

  // Check that all parameters were set correctly
  this->CheckParameters();

  // Generate the dynamic image pyramids
  for ( unsigned int i = 0; i < this->numberOfDynamicImages; ++i )
  {
    auto curDynPyramid = std::make_shared<ImagePyramid<VoxelType> >();
    curDynPyramid->GenerateLevels( this->allDynamicImages[i], this->numberOfLevels, this->numberOfLevelsToPerform );
    allDynamicPyramids.push_back( curDynPyramid );
  }

  // Generate the reference state image pyramid
  this->referenceStatePyramid = std::make_shared<ImagePyramid<VoxelType> >();
  this->referenceStatePyramid->GenerateLevels( this->referenceStateImage,
                                               this->numberOfLevels,
                                               this->numberOfLevelsToPerform,
                                               2 );
  

  // Initialise the transformation
  if (this->transformationType == STANDARD_B_SPLINE)
  {
    auto bsplTrafo = std::make_shared<BSplineTransformation>( this->defSpaceImage,
      this->numberOfLevelsToPerform,
      this->bSplineCPGSpacing );
    bsplTrafo->SetBendingEnergyWeight( this->bSplineBendingEnergyWeight );
    bsplTrafo->SetLinearEnergyWeight( this->bSplineLinearEnergyWeight );
    this->transform = bsplTrafo;
  }
  else if (SLIDING_B_SPLINE == this->transformationType)
  {
    // We need two transforms for setting up the sliding transform, these do not have to be the same
    auto bsplTrafoA = std::make_shared<BSplineTransformation>( this->defSpaceImage,
      this->numberOfLevelsToPerform,
      this->bSplineCPGSpacing );
    bsplTrafoA->SetBendingEnergyWeight( this->bSplineBendingEnergyWeight );
    bsplTrafoA->SetLinearEnergyWeight( this->bSplineLinearEnergyWeight );
    auto bsplTrafoB = bsplTrafoA->DeepCopy();

    // Use the two transformations above to generate a sliding transformation
    auto slidingBSplineTrafo = std::make_shared<SlidingTransformation>( bsplTrafoA, bsplTrafoB, this->defSpaceImage, numberOfLevels, numberOfLevelsToPerform );
    slidingBSplineTrafo->SetGapOverlapConstraintWeight( this->slidingGapOverlapConstraintWeight );
    slidingBSplineTrafo->SetSignedDistanceMapImage( this->slidingSignedDistanceMapImage );

    this->transform = slidingBSplineTrafo;
  }
  
  // Initialise the correspondence model
  this->correspondenceModel = std::make_shared<CorrespondenceModel>( this->numberOfSurrogateSignals, this->transform );

  // Check if an input correspondence model was provided
  if (this->inputRCMImages.size() > 0)
  {
    std::cout << "Input correspondence model given. Findin appropriate multi-resolution level..." << std::endl;
    
    // Check how many voxels are in the given images
    unsigned int totalNumOfInputRCMParameters = 0;
    
    for (size_t iRCMImg = 0; iRCMImg < this->inputRCMImages.size(); ++iRCMImg)
    {
      totalNumOfInputRCMParameters += this->inputRCMImages[iRCMImg]->nvox;
    }
    std::cout << "Total number of input correspondence model parameters: " << totalNumOfInputRCMParameters << std::endl;

    for (this->currentLevel = 0; this->currentLevel < numberOfLevelsToPerform; ++this->currentLevel)
    {
      this->correspondenceModel->InitialiseLevel( this->currentLevel );
      /// Todo: Remove transform. Does not need to be a member variable!
      this->transform->InitialiseLevel( this->currentLevel ); 
      unsigned int numberOfCorrespondenceModelParameters = this->correspondenceModel->GetNumberOfParameters();
      
      if (numberOfCorrespondenceModelParameters != totalNumOfInputRCMParameters)
      {
        std::cout << "Mismatch of parameter numer on level: " << this->currentLevel
          << " input: " << totalNumOfInputRCMParameters
          << " vs. expected: " << numberOfCorrespondenceModelParameters << std::endl;
        continue;
      }
      else
      {
        std::cout << "Match of parameter numer on level: " << this->currentLevel << std::endl;
        std::cout << "Processing will commence on this level: " << this->currentLevel << std::endl;
        break;
      }
    }

    this->correspondenceModel->SetParameters( this->inputRCMImages );
  }

  // Initialise the image similarity measure
  // Currently only SSD supported 
  this->similarityMeasure = std::make_shared<SSDImageSimilarity>();

  // Initialise the image acquisition
  switch (this->dynamicImageDataType)
  {
    case SAME_RES_AS_STATIC:
    {
      auto noImageAcquisition = std::make_shared<NoImageAcquisition>();
      this->imageAcquisition = noImageAcquisition;
      break;
    }
    case LOWER_RES_THAN_STATIC:
    {
      auto lowResImageAcquisition = std::make_shared<LowResolutionImageAcquisition>();
      this->imageAcquisition = lowResImageAcquisition;
      break;
    }
    default:
      // Defaults to no-image-acquisition simulation
      auto noImageAcquisition = std::make_shared<NoImageAcquisition>();
      this->imageAcquisition = noImageAcquisition;
      break;
  }

  this->initialised = true;
  }




//--------------------------
// Supremo::CheckParameters
//--------------------------
void Supremo::CheckParameters()
{
  // Check the reference image
  if (this->referenceStateImage == nullptr)
  {
    char msg[200];
    sprintf_s( msg, "Reference state image was not set." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Check all dynamic images
  for ( unsigned int i = 0; i < this->numberOfDynamicImages; ++i )
  {
    if ( this->allDynamicImages[i] == nullptr )
    {
      char msg[200];
      sprintf_s( msg, "Dynamic image %i was not set.", i );
      supremo_print_error( msg );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }

  // Check if the image that defines the space was set, otherwise use the reference state image
  if (this->defSpaceImage == nullptr)
  {
    this->defSpaceImage = this->referenceStateImage;
  }

  // number of levels and number of levels to perform
  if ( this->numberOfLevels <= 0 )
  {
    char msg[200];
    sprintf_s( msg, "Number of levels was set to %i, but has to be larger than zero.", this->numberOfLevels );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  if ( this->numberOfLevelsToPerform <= 0 )
  {
    char msg[200];
    sprintf_s( msg, "Number of levels to perform was set to %i, but has to be larger than zero.", this->numberOfLevelsToPerform );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Check the b-spline control-point spacing
  // Set the spacing along y and z if undefined. Their values are set to match
  // the spacing along the x axis
  if (this->bSplineCPGSpacing[0] != this->bSplineCPGSpacing[0])
  {
    char msg[200];
    sprintf_s( msg, "B-spline control point spacing was set to undefined value." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  if (this->bSplineCPGSpacing[1] != this->bSplineCPGSpacing[1]) this->bSplineCPGSpacing[1] = this->bSplineCPGSpacing[0];
  if (this->bSplineCPGSpacing[2] != this->bSplineCPGSpacing[2]) this->bSplineCPGSpacing[2] = this->bSplineCPGSpacing[0];

  // Check that the number of surrogate signals is geq 1
  if ( this->numberOfSurrogateSignals <= 1 )
  {
    char msg[200];
    sprintf_s( msg, "Number of surrogate signals must not be lower than 1." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  //check if maxSwitchIterationNumber has been set (will be -1 if not)
  //if not set to 10 if motion compensated reconstruction is being used, to 1 if not
  if (this->maxSwitchIterationNumber < 0)
  {
    if (this->mocoReconstructor == nullptr)
      this->maxSwitchIterationNumber = 1;
    else
      this->maxSwitchIterationNumber = 10;
  }
}




//----------------------------------------
// Supremo::DisplayCurrentLevelParameters
//----------------------------------------
void Supremo::DisplayCurrentLevelParameters(unsigned int levelIn)
{
  printf( "Current level: %i / %i\n", levelIn+1, this->numberOfLevelsToPerform );
  
  // Display the reference state image parameters
  printf( "Reference state image\n");
  printf( "\t* image dimension: %i x %i x %i\n", 
          this->currentReferenceStateImage->nx,
          this->currentReferenceStateImage->ny,
          this->currentReferenceStateImage->nz );
  printf( "\t* image spacing: %g x %g x %g mm\n\n", 
          this->currentReferenceStateImage->dx,
          this->currentReferenceStateImage->dy,
          this->currentReferenceStateImage->dz );
  
  // Display the transoformation parameters
  this->transform->DisplayTransformationParameters();

  // Display the dynamic image parameters
  printf( "1st Dynamic image\n");
  printf( "\t* image dimension: %i x %i x %i\n",
          this->allDynamicPyramids[0]->GetLevel( levelIn )->nx,
          this->allDynamicPyramids[0]->GetLevel( levelIn )->ny,
          this->allDynamicPyramids[0]->GetLevel( levelIn )->nz );
  printf( "\t* image spacing: %g x %g x %g mm\n",
          this->allDynamicPyramids[0]->GetLevel( levelIn )->dx,
          this->allDynamicPyramids[0]->GetLevel( levelIn )->dy,
          this->allDynamicPyramids[0]->GetLevel( levelIn )->dz );

#ifdef _DEBUG


  for (unsigned int n = 0; n < this->numberOfDynamicImages; ++n)
  {
    
    printf( "Dynamic image: %i\n", n );

    printf( "\t* image dimension: %i x %i x %i\n",
            this->allDynamicPyramids[n]->GetLevel( levelIn )->nx,
            this->allDynamicPyramids[n]->GetLevel( levelIn )->ny,
            this->allDynamicPyramids[n]->GetLevel( levelIn )->nz );

    printf( "\t* image spacing: %g x %g x %g\n",
            this->allDynamicPyramids[n]->GetLevel( levelIn )->dx,
            this->allDynamicPyramids[n]->GetLevel( levelIn )->dy,
            this->allDynamicPyramids[n]->GetLevel( levelIn )->dz );

    if (this->allDynamicPyramids[n]->GetLevel(levelIn)->sform_code>0)
      reg_mat44_disp( &(this->allDynamicPyramids[n]->GetLevel( levelIn )->sto_xyz), (char *)"\t* sform matrix" );
    else reg_mat44_disp( &(this->allDynamicPyramids[n]->GetLevel( levelIn )->qto_xyz), (char *)"\t* qform matrix" );
  }
#endif

  return;
}





void Supremo::SaveInterMCRToImage( nifti_image* mcrImage )
{  
  // write the reconstructed image to the defined destination
  if (!this->outputInterMCRFolder.empty())
  {
    std::ostringstream osOutFileName;
    osOutFileName << this->outputInterMCRFolder << "/MCR_" << this->levelID << ".nii.gz";
    osOutFileName.str().c_str();
    reg_io_WriteImageFile( mcrImage, osOutFileName.str().c_str() );
  }
}

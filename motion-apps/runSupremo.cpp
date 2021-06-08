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



#include <memory>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>  
#include "nifti1_io.h"
#include "_reg_ReadWriteImage.h"
#include "Supremo.h"
#include "GitSHA1.h"

// Supremo specific includes
#include "ImagePyramid.h"
#include "CommandLineParser.h"
#include "MoCoRecon.h"
#include "MoCoReconWeightedAverage.h"
#include "MoCoReconSuperResolutionIBP.h"




/** 
 *  Main executable that instantiates Supremo
 */
int main( int argc, char *argv[] )
{
  // First things first: Which version is actually running:
  std::cout << "Running Supremo version (git commit): " << g_GIT_SHA1 << std::endl;
  
  
  // Specify a map that lists all known parameters
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-refState"] = { 1, true, "Reference state image" };
  commandLineOptions["-dynamic"]  = { 2, true, "Dynamic image data" };
  commandLineOptions["-surr"]     = { 2, true, "Surrogate data" };

  commandLineOptions["-dType"]        = { 1, false, "Dynamic image type. Determines the image acuqisiton simulation" };
  commandLineOptions["-defSpace"]     = { 1, false, "Defines the space of the dynamic image" };
  commandLineOptions["-mcrType"]      = { 1, false, "Type of motion compensated imgae reconstruction" };
  commandLineOptions["-maxMCRIt"]     = { 1, false, "Maximum number of MCR iterations" };
  commandLineOptions["-outRCM"]       = { 1, false, "Output file name for respiratory correspondence model" };
  commandLineOptions["-outSimDyn"]    = { 1, false, "Output file directory for the simulated dynamic images" };
  commandLineOptions["-outDVFs"]      = { 1, false, "Output file directory for the correspondence model-generated DVFS. The DVFs will have the size of the dynamic images." };
  commandLineOptions["-outMCR"]       = { 1, false, "Output file name for the motion compensated reconstructed image" };
  commandLineOptions["-outInterMCR"]  = { 1, false, "Output directory name for the intermediate motion compensated reconstructed images" };
  commandLineOptions["-outInterGrad"] = { 1, false, "Output directory name for the intermediate gradients of the objective function" };
  commandLineOptions["-inRCM"]        = { 1, false, "Input name for a respiratory correspondence model. Comma-separated list expected" };
  commandLineOptions["-transType"]    = { 1, false, "Transformation type" };
  commandLineOptions["-distMap"]      = { 1, false, "Input distance map" };
  commandLineOptions["-sx"]           = { 1, false, "B-spline control point spacing in x-direciton" };
  commandLineOptions["-sy"]           = { 1, false, "B-spline control point spacing in y-direciton" };
  commandLineOptions["-sz"]           = { 1, false, "B-spline control point spacing in z-direciton" };
  commandLineOptions["-be"]           = { 1, false, "Bending energy weight" };
  commandLineOptions["-le"]           = { 1, false, "Linear energy weight" };
  commandLineOptions["-go"]           = { 1, false, "Gap/overlap constraint weight when using sliding transformation" };
  commandLineOptions["-maxSwitchIt"]  = { 1, false, "Maximum nuimber of switch iterations" };
  commandLineOptions["-ln"]           = { 1, false, "Number of pyramid levels" };
  commandLineOptions["-lp"]           = { 1, false, "Number of pyramid levels to be performed" };
  commandLineOptions["-maxFitIt"]     = { 1, false, "Maximum number of fit iterations" };
  commandLineOptions["-h"]            = { 1, false, "Print help message." };

  // Parse the command line
  std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>( argc, argv, commandLineOptions );
  std::cout << parser->getCommandLine() << std::endl;

  // Check if help was required 
  if (parser->cmdOptionExists( "-h" ) || parser->getAllReqreuiredParametersSet())
  {
      // TODO: print help and exit
      std::cout << commandLineOptions["-refState"].description << std::endl;
      return EXIT_FAILURE;
  }

  //--------------------------------
  // Load the reference-state image 
  //--------------------------------
  /// \todo: invalid file name causes memory access violation
  std::string referenceStateImageFileName = parser->getCmdOptionAsString( "-refState" );
  nifti_image* referenceStateImage = nifti_image_read( referenceStateImageFileName.c_str(), true );

  if (referenceStateImage == nullptr)
  {
    char msg[200];
    sprintf_s( msg, "Could not read reference state image: %s", parser->getCmdOptionAsString( "-refState" ).c_str() );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  //-------------------------
  // Load the dynamic images 
  //-------------------------

  int numberOfDynamicImages = parser->getCmdOptionAsInt( "-dynamic", 0 );
  std::string dynamicImageFileName = parser->getCmdOptionAsString( "-dynamic", 1 );
  nifti_image** dynamicImages = nullptr; // (nifti_image **)malloc( numberOfDynamicImages * sizeof( nifti_image * ) );
  std::ifstream dynamicNamesFile( dynamicImageFileName.c_str(), std::ifstream::in );

  if (!dynamicNamesFile.is_open())
  {
    char msg[200];
    sprintf_s( msg, "Cannot open the dynamic image names file %s", dynamicImageFileName.c_str() );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Read the complete file, check if correct number of files was specified in the file, then read
  {
    std::vector<std::string> allDynamicImageNames;
    std::string curDynamicImageName;

    // read values until no more a found in the file
    while (dynamicNamesFile >> curDynamicImageName)
    {
      allDynamicImageNames.push_back( curDynamicImageName );
    }

    // Check that the correct number of dynamic images was provided
    if (allDynamicImageNames.size() != numberOfDynamicImages)
    {
      supremo_print_error( "Number of dynamic images not as expected." );
      supremo_exit( 1, __FILE__, __LINE__ );
    }

    // Now read the images
    dynamicImages = (nifti_image **)malloc( numberOfDynamicImages * sizeof( nifti_image * ) );

    for (int d = 0; d < numberOfDynamicImages; ++d)
    {
      dynamicImages[d] = nifti_image_read( allDynamicImageNames[d].c_str(), true );
      
      // Check if the current image was loaded properly
      if (dynamicImages[d] == nullptr)
      {
        char msg[200];
        sprintf_s( msg, "Unable to open dynamic image %i: %s", d, allDynamicImageNames[d].c_str() );
        supremo_print_error( msg );
        supremo_exit( 1, __FILE__, __LINE__ );
      }
    }
  }
  dynamicNamesFile.close();


  //-------------------------
  // Load the surrogate data
  //-------------------------

  int numberOfSurrogateSignals = parser->getCmdOptionAsInt( "-surr", 0 );
  std::string surrogateFileName = parser->getCmdOptionAsString( "-surr", 1 );

  // Open the file  
  std::ifstream surrSignalFile( surrogateFileName.c_str(), std::ifstream::in );

  // Check file was opened correctly
  if (!surrSignalFile.is_open())
  {
    char msg[200];
    sprintf_s( msg, "Surrogate signal file could not be opened: %s", surrogateFileName.c_str() );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // generate a variable for the surrogate signals in this scope
  float* surrogateSignals = nullptr;

  // Read the complete surrogate file first and proceed only if the size was as expected
  {
    float curSurrVal;
    std::vector<float> rawSurrValues;

    // read values until no more a found in the file
    while (surrSignalFile >> curSurrVal)
    {
      rawSurrValues.push_back( curSurrVal );
    }

    // Check that the correct number of surrogate signals was provided
    if (rawSurrValues.size() != numberOfSurrogateSignals * numberOfDynamicImages)
    {
      supremo_print_error( "Number of surrogate signals not as expected." );
      supremo_exit( 1, __FILE__, __LINE__ );
    }

    // allocate float array of correct size and copy over values from the vector
    // we could just take a pointer like ( &rawSurrValues[0] ), but dynamic memory
    // management of the vector might invalidate that pointer. Hence explicit.
    surrogateSignals = new float[numberOfSurrogateSignals * numberOfDynamicImages];

    for (int i = 0; i < numberOfDynamicImages*numberOfSurrogateSignals; ++i)
    {
      surrogateSignals[i] = rawSurrValues[i];
    }
  }
  surrSignalFile.close();

  //------------------------------------
  // Get remaining command line options
  //------------------------------------

  // -defSpace
  nifti_image* defSpaceImage = nullptr;

  if (parser->cmdOptionExists( "-defSpace" ))
  {
    defSpaceImage = reg_io_ReadImageFile( parser->getCmdOptionAsString( "-defSpace" ).c_str() );

    if (defSpaceImage == nullptr)
    {
      char msg[200];
      sprintf_s( msg, "Could not read space definition image: %s", parser->getCmdOptionAsString( "-defSpace" ).c_str() );
      supremo_print_error( msg );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }

  // -dType
  // Default to "same resolution" option
  t_dynamicData dynImageType = SAME_RES_AS_STATIC;
  
  if (parser->cmdOptionExists( "-dType" ))
  {
    int iTmpDynType = parser->getCmdOptionAsInt( "-dType" );
    switch (iTmpDynType)
    {
      case SAME_RES_AS_STATIC:
        dynImageType = SAME_RES_AS_STATIC;
        break;
      case LOWER_RES_THAN_STATIC:
        dynImageType = LOWER_RES_THAN_STATIC;
        break;
      default:
        char msg[200];
        sprintf_s( msg, "Dynamic image type unknown: %i", iTmpDynType );
        supremo_print_error( msg );
        supremo_exit( 1, __FILE__, __LINE__ );
    }
  }

  // -mcrType
  // default to "no reconstruction" option
  // Build the reconstruction object here and pass it to supremo
  std::shared_ptr<MoCoRecon> motionCompensatedImageRecon = nullptr;

  //// Motion compensated image reconstruction 
  if (parser->cmdOptionExists( "-mcrType" ))
  {
    int iTmpMCRType = parser->getCmdOptionAsInt( "-mcrType" );
    
    // -maxMCRIt
    // default to 5
    unsigned int maxMCRIt = 5;
    if (parser->cmdOptionExists( "-maxMCRIt" ))
    {
      maxMCRIt = static_cast<unsigned int>(parser->getCmdOptionAsInt( "-maxMCRIt" ));
    }

    switch (iTmpMCRType)
    {
      case NO_RECONSTRUCTION:
      {
        // leave the pointer to be a null pointer if no image reconstruction is required
        motionCompensatedImageRecon = nullptr;
        break;
      }
      case WEIGHTED_AVERAGING:
      {
        auto avgWeightedReconstructor = std::make_shared<MoCoReconWeightedAverage>();
        motionCompensatedImageRecon = avgWeightedReconstructor;
        break;
      }
      case SUPER_RESOLUTION_RESTART:
      {
        auto superResReconstructor = std::make_shared<MoCoReconSuperResolutionIBP>( maxMCRIt, false );
        motionCompensatedImageRecon = superResReconstructor;
        break; 
      }
      case SUPER_RESOLUTION_UPDATE:
      {
        auto superResReconstructor = std::make_shared<MoCoReconSuperResolutionIBP>( maxMCRIt, true );
        motionCompensatedImageRecon = superResReconstructor;
        break;
      }
      default:
      {
        char msg[200];
        sprintf_s( msg, "Motion compensated image reconstruction type unknown: %i", iTmpMCRType );
        supremo_print_error( msg );
        supremo_exit( 1, __FILE__, __LINE__ );
      }
    }
  }

  // -outInterMCR
  // default to empty string
  std::string outputInterMCRFolderName = "";
  if (parser->cmdOptionExists( "-outInterMCR" ))
  {
    outputInterMCRFolderName = parser->getCmdOptionAsString( "-outInterMCR" );
  }
  
  // -outInterGrad
  // default to empty string
  std::string outputInterGradientFolderName = "";
  if (parser->cmdOptionExists( "-outInterGrad" ))
  {
    outputInterGradientFolderName = parser->getCmdOptionAsString( "-outInterGrad" );
  }


  // -transType
  // default to standard B-spline
  t_transformationType transformType = STANDARD_B_SPLINE;
  if (parser->cmdOptionExists( "-transType" ))
  {
    int tmpTransType = parser->getCmdOptionAsInt( "-transType" );

    switch (tmpTransType)
    {
    case STANDARD_B_SPLINE:
      transformType = STANDARD_B_SPLINE;
      break;
    case SLIDING_B_SPLINE:
      transformType = SLIDING_B_SPLINE;
      break;
    default:
      char msg[200];
      sprintf_s( msg, "Transformation type unknown: %i", tmpTransType );
      supremo_print_error( msg );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }

  // -distMap
  // Singed distance map required for sliding registartion
  nifti_image* slidingSingedDistanceMap = nullptr;
  if (parser->cmdOptionExists( "-distMap" ))
  {
    slidingSingedDistanceMap = nifti_image_read( parser->getCmdOptionAsString( "-distMap" ).c_str(), true );
    
    // Check that the image was read successfully 
    if (slidingSingedDistanceMap == nullptr)
    {
      char msg[200];
      sprintf_s( msg, "Could not read signed distance map required for sliding transformation: %s", parser->getCmdOptionAsString( "-distMap" ).c_str() );
      supremo_print_error( msg );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }
  
  // -sx, -sy, -sz
  // default to 5 voxels, i.e. sx<0
  float sx=-5, sy, sz;
  if (parser->cmdOptionExists( "-sx" ))  sx = parser->getCmdOptionAsFloat( "-sx" );
  if (parser->cmdOptionExists( "-sy" ))  sy = parser->getCmdOptionAsFloat( "-sy" );
  else sy = sx;
  if (parser->cmdOptionExists( "-sz" ))  sz = parser->getCmdOptionAsFloat( "-sz" );
  else sz = sx;

  // -be
  // default to 0
  float bSplineBendingEnergyWeight = 0;
  if (parser->cmdOptionExists( "-be" )) bSplineBendingEnergyWeight = parser->getCmdOptionAsFloat( "-be" );

  // -le
  // default to 0
  float bSplineLinearEnergyWeight = 0;
  if (parser->cmdOptionExists( "-le" )) bSplineLinearEnergyWeight = parser->getCmdOptionAsFloat( "-le" );

  // -go
  // default to 0
  float slidingGapOverlapConstraintWeight = 0;
  if (parser->cmdOptionExists( "-go" )) slidingGapOverlapConstraintWeight = parser->getCmdOptionAsFloat( "-go" );

  // -ln
  //defaults to 3
  unsigned int numberOfLevels = 3;
  if (parser->cmdOptionExists( "-ln" )) numberOfLevels = parser->getCmdOptionAsInt( "-ln" );

  // -lp
  //defaults to number of levels
  unsigned int numberOfLevelsToPerform = numberOfLevels;
  if (parser->cmdOptionExists( "-lp" )) numberOfLevelsToPerform = parser->getCmdOptionAsInt( "-lp" );

  // -maxFitIt
  // defaults to 300
  unsigned int maxModelFittingIterations = 50;
  if (parser->cmdOptionExists( "-maxFitIt" )) maxModelFittingIterations = static_cast<unsigned int>(parser->getCmdOptionAsInt( "-maxFitIt" ));

  // -maxSwitchIt
  unsigned int maxSwitchIterations = 10;
  if (parser->cmdOptionExists( "-maxSwitchIt" )) maxSwitchIterations = static_cast<unsigned int>(parser->getCmdOptionAsInt( "-maxSwitchIt" ));


  // -inRCM
  // Read the input respiratory correspondence model image(s)
  std::vector<nifti_image*> inputRCMImages;
  if (parser->cmdOptionExists( "-inRCM" ))
  {
    std::vector<std::string> rcmImageFileNames = splitStringbyDelimiter( parser->getCmdOptionAsString( "-inRCM" ), "," );
    
    for (size_t nImg = 0; nImg < rcmImageFileNames.size(); ++nImg)
    {
      nifti_image* curRCMInImg = nifti_image_read( rcmImageFileNames[nImg].c_str(), true );
      if (nullptr != curRCMInImg)
      {
        inputRCMImages.push_back( curRCMInImg );
      }
      else
      {
        char msg[200];
        sprintf_s( msg, "Could not read input RCM image: %s", rcmImageFileNames[nImg].c_str() );
        supremo_print_error( msg );
        supremo_exit( 1, __FILE__, __LINE__ );
      }
    }
  }


  //------------------------------------------------
  // Feed Supremo with the input data and options
  //------------------------------------------------
  std::shared_ptr<Supremo> supremo = std::make_shared<Supremo>();
  supremo->SetDynamicImages( dynamicImages, numberOfDynamicImages );
  supremo->SetReferenceStateImage( referenceStateImage );
  supremo->SetSurrogateSignals( surrogateSignals, numberOfSurrogateSignals );
  supremo->SetDefSpaceImage( defSpaceImage ); 
  supremo->SetDynamicImageDataType( dynImageType );
  supremo->SetInputRCMImages( inputRCMImages );
  supremo->SetMotionCompensatedReconstruction( motionCompensatedImageRecon );
  supremo->SetInterMCROutputFolder( outputInterMCRFolderName );
  supremo->SetInterGradOutputFolder( outputInterGradientFolderName );
  supremo->SetBSplineCPGSpacing( sx, sy, sz );
  supremo->SetBSplineBendingEnergy( bSplineBendingEnergyWeight );
  supremo->SetBSplineLinearEnergy( bSplineLinearEnergyWeight );
  supremo->SetTransformationType( transformType );
  supremo->SetSlidingTrafoDistanceMap( slidingSingedDistanceMap );
  supremo->SetNumberOfPyramidLevels( numberOfLevels );
  supremo->SetNumberOfPyramidLevelsToPerform( numberOfLevelsToPerform );
  supremo->SetMaxModelFittingIterationNumber( maxModelFittingIterations );
  supremo->SetMaxSwitchIterationNumber( maxSwitchIterations );
  supremo->FitMotionModelAndReconstruct();
  

  
  //--------------
  // Save results
  //--------------

  // Respiratory correspondence model images will always be saved. Define a name if this was not done. 
  std::string outputRCMFileName = "outputRCM.nii.gz";
  
  if (parser->cmdOptionExists( "-outRCM" ))
  {
    outputRCMFileName = parser->getCmdOptionAsString( "-outRCM" );
  }

  // Remove anything after .nii (this includes .nii.gz) from the required file name such that we can append the index of the images later
  outputRCMFileName.replace( outputRCMFileName.find( ".nii" ), outputRCMFileName.length(), "" );

  // Get the images and save each one individually
  std::vector<nifti_image*> respCorrModelImgs = supremo->GetCorrespondenceModelAsImage();
  
  for (size_t i = 0; i < respCorrModelImgs.size(); ++i)
  {
    std::ostringstream osOutFileName;
    osOutFileName << outputRCMFileName << "_t" << std::setfill( '0' ) << std::setw( 2 ) << i << ".nii.gz";

    reg_io_WriteImageFile( respCorrModelImgs[i], osOutFileName.str().c_str() );
    nifti_image_free( respCorrModelImgs[i] );
    respCorrModelImgs[i] = nullptr;
  }
    
  // Save the motion-compensated recontructed image if required
  if (parser->cmdOptionExists( "-outMCR" ))
  {
    outputRCMFileName = parser->getCmdOptionAsString( "-outMCR" );

    nifti_image* mocoReconImage = supremo->GetMotionCompensatedReconstructedImage();
    reg_io_WriteImageFile( mocoReconImage, parser->getCmdOptionAsString( "-outMCR" ).c_str() );
  }
  
  // Save the simulated dynamic images
  if (parser->cmdOptionExists( "-outSimDyn" ))
  {
    std::string outSimDynFolder = parser->getCmdOptionAsString( "-outSimDyn" );
    
    std::vector<nifti_image*> simulatedDynImgs = supremo->SimulateDynamicImages();

    for (unsigned int i = 0; i < simulatedDynImgs.size(); ++i)
    {
      char fOutName[200];
      sprintf( fOutName, "%ssimDyn%04i.nii.gz",  outSimDynFolder.c_str(), i  );
      reg_io_WriteImageFile( simulatedDynImgs[i], fOutName );

      // Directly clean up the simulated dynamic images
      nifti_image_free( simulatedDynImgs[i] );
      simulatedDynImgs[i] = nullptr;
    }
  }
  
  // Save the deformation vector fields for the individual time points
  if (parser->cmdOptionExists( "-outDVFs" ))
  {
    std::string outDVFFolder = parser->getCmdOptionAsString( "-outDVFs" );

    std::vector<nifti_image*> generatedDVFs = supremo->GenerateDVFsFromCorrespondenceModel();

    for (unsigned int i = 0; i < generatedDVFs.size(); ++i)
    {
      char fOutName[200];
      sprintf( fOutName, "%sDVF%04i.nii.gz", outDVFFolder.c_str(), i );
      reg_io_WriteImageFile( generatedDVFs[i], fOutName );

      // Free up this image
      nifti_image_free( generatedDVFs[i] );
      generatedDVFs[i] = nullptr;
    }
  }



  //if (parser->cmdOptionExists( "-outSliDistMap" ))
  //{
  //  std::string outSlidingDistMapOutputFile = parser->getCmdOptionAsString( "-outSliDistMap" );
  //  // ToDo: Save final map used for sliding transformation
  //}

  //if (parser->cmdOptionExists( "-outSliShapeModelParams" ))
  //{
  //  std::string outSlidingDistMapOutputFile = parser->getCmdOptionAsString( "-outSliShapeModelParams" );
  //  // ToDo: Save final motion mask parametersdistance map used for sliding transformation
  //}


  


  //----------
  // Clean up
  //----------

  // Delete the Supremo object
  supremo.reset();

  // Delete the command line parser
  parser.reset();

  // Clean up dynamic images if they were loaded by Supremo (as opposed to loading them from outside and cleaning up there)
  if (dynamicImages != nullptr)
  {
    // clean up each image 
    for (int i = 0; i < numberOfDynamicImages; i++)
    {
      if (dynamicImages[i] != nullptr) nifti_image_free( dynamicImages[i] );
    }
    // then free the memory
    free( dynamicImages );
  }


  // Clean up the reference state image.
  if (referenceStateImage != nullptr)
  {
    nifti_image_free( referenceStateImage );
    referenceStateImage = nullptr;
  }

  // Clean up the surrogate data
  if (surrogateSignals != nullptr)
  {
    delete[] surrogateSignals;
  }

  return EXIT_SUCCESS;
}

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
int main(int argc, char* argv[])
{
  // First things first: Which version is actually running:
  std::cout << "Running Supremo version (git commit): " << g_GIT_SHA1 << std::endl;


  // Specify a map that lists all known parameters
  std::map<std::string, CommandLineOption> commandLineOptions;
  // Required
  commandLineOptions["-refState"] = { 1, true, "Reference state image (nifti format)", "<filename>" };
  commandLineOptions["-dynamic"] = { 2, true, "Dynamic image data (Number of dynamic images and name of a text "
                                               "file containing the list of dynamic images)", "<int> <filename>" };
  commandLineOptions["-surr"] = { 2, true, "Surrogate data (Number of surrogate signals and name of a text "
                                               "file containing the values of the surrogate signal(s).These should "
                                               "be ordered with 1st value of each signal, then 2nd value of each "
                                               "signal, ..., then last value of each signal)", "<int> <filename>" };
  // Dynamic image data
  commandLineOptions["-dType"] = { 1, false, "Specify type of dynamic images, where: \n"
                                                    "[0] = Full resolution images(default) \n"
                                                    "Can be full or partial image data(slices, slabs) that has the same "
                                                    "resolution as the reference state image.The nifti file headers should "
                                                    "provide the imaging geometry.\n"
                                                    "[1] = Low resolution images\n"
                                                    "Can be full or partial image data(slices, slabs) that has a lower "
                                                    "resolution than the reference state image.Acquisition of low res data will be "
                                                    "simulated along any dimensions where dynamic voxel size > reference state voxel "
                                                    "size. The nifti file headers should provide the imaging geometry.", "<int>" };
  commandLineOptions["-defSpace"] = { 1, false, "Filename of the image used to define the space of the deformed images. "
                                                    "This image is used to define the extent and resolution(if specified in "
                                                    "voxels) of the modeland transform CPGs.If a defSpace image is not "
                                                    "specified the reference state image will be used.", "<filename>" };
  // Motion-compensated image reconstruction
  commandLineOptions["-mcrType"] = { 1, false, "Specify type of motion compensated image reconstruction, where:\n"
                                                    "[0] = No motion compensated image reconstruction(default)\n"
                                                    "[1] = Weighted average of deformed dynamic images\n"
                                                    "[2] = Super resoltuion(iterative back - projection) - restart recon\n"
                                                    "[3] = Super resoltuion(iterative back - projection) - update recon", "<int>" };
  commandLineOptions["-maxMCRIt"] = { 1, false, "Maximum number of iterations to use with iterative reconstruction methods [5]\n"
                                                    "When motion compensated image reconstruction is performed the reference state image "
                                                    "provided as input is used to define the space of the reconstructed image.", "<int>" };
  // Out-/input options
  commandLineOptions["-outRCM"] = { 1, false, "Filename for saving the respiratory correspondence model [outputRCM.nii.gz]", "<filename>" };
  commandLineOptions["-outMCR"] = { 1, false, "The final motion compensated image reconstruction will be saved using the "
                                                    "filename provided", "<filename>" };
  commandLineOptions["-outSimDyn"] = { 1, false, "The final simualted dynamic image data will be saved in the folder specified", "<folder>" };
  commandLineOptions["-outDVFs"] = { 1, false, "The final correspondence model-generated DVFs will be saved in the folder "
                                                    "specified.The DVFs will have the size of the dynamic images.", "<folder>" };
  commandLineOptions["-outInterMCR"] = { 1, false, "The intermediate MCRs will be saved to the folder specified", "<folder>" };
  commandLineOptions["-outInterGrad"] = { 1, false, "The intermediate objective function gradients will be saved to the folder specified", "<folder>" };
  commandLineOptions["-inRCM"] = { 1, false, "Input file name(s) for a respiratory correspondence model (comma-separated list, one file per sliding "
                                                    "region required if used).", "<fName1>,<fName2>" };
  // Transformation options
  commandLineOptions["-transType"] = { 1, false, "Transformation type [0]\n"
                                                    "[0] = B-spline transformation\n"
                                                    "[1] = Sliding B-spline transformation(requires signed distance map)", "<int>" };
  commandLineOptions["-sx"] = { 1, false, "Final grid spacing along the x axis in mm (in voxel if negative value) [5 voxels]", "<float>" };
  commandLineOptions["-sy"] = { 1, false, "Final grid spacing along the y axis in mm (in voxel if negative value) [sx value]", "<float>" };
  commandLineOptions["-sz"] = { 1, false, "Final grid spacing along the z axis in mm (in voxel if negative value) [sx value]", "<float>" };
  commandLineOptions["-be"] = { 1, false, "Weight of the bending energy penalty term [0.0]", "<float>" };
  commandLineOptions["-le"] = { 1, false, "Weight of the first order penalty term (symmetric and anti-symmetric part of the Jacobian) [0.0]", "<float>" };
  commandLineOptions["-distMap"] = { 1, false, "Signed distance map defining boundary of sliding regions at the zero-crossing", "<filename>" };
  commandLineOptions["-go"] = { 1, false, "Weight of the gap-overlap penalty term when using sliding transformation", "<float>" };
  // Optimisation options
  commandLineOptions["-maxSwitchIt"] = { 1, false, "Maximum number of times to iterate between motion compensate image reconstruction and fitting "
                                                    "the respiratory correspondence model[10]", "<int>" };
  commandLineOptions["-ln"] = { 1, false, "Number of level imagepyramid levels to generate [3]", "<int>" };
  commandLineOptions["-lp"] = { 1, false, "Only perform processing on the first lp levels [ln]", "<int>" };
  commandLineOptions["-maxFitIt"] = { 1, false, "Maximum number of respiratory correspondence model fitting iterations [300]", "<int>" };
  // Help/information
  commandLineOptions["-h"] = { 0, false, "Print help message and exit.", "" };


  // Parse the command line
  std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>(argc, argv, commandLineOptions);
  std::cout << parser->getCommandLine() << std::endl;

  // Check if help is required 
  if (parser->cmdOptionExists("-h") || !parser->getAllReqreuiredParametersSet())
  {
    // Required: 
    std::vector<std::string> opts;
    opts.push_back("-refState");
    opts.push_back("-dynamic");
    opts.push_back("-surr");
    printFormattedCommandLineOptions(commandLineOptions, "Required options", 100, opts);

    // Dynamic image data options
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-dType");
    opts.push_back("-defSpace");
    printFormattedCommandLineOptions(commandLineOptions, "Dynamic image data options", 100, opts);

    // Motion-compensated image reconstruction
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-mcrType");
    opts.push_back("-maxMCRIt");
    printFormattedCommandLineOptions(commandLineOptions, "Motion-compensated image reconstruction", 100, opts);

    // Out-/input options
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-outRCM");
    opts.push_back("-outMCR");
    opts.push_back("-outSimDyn");
    opts.push_back("-outDVFs");
    opts.push_back("-outInterMCR");
    opts.push_back("-outInterGrad");
    opts.push_back("-inRCM");
    printFormattedCommandLineOptions(commandLineOptions, "Out-/input options", 100, opts);

    // Transformation options
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-transType");
    opts.push_back("-sx");
    opts.push_back("-sy");
    opts.push_back("-sz");
    opts.push_back("-be");
    opts.push_back("-le");
    opts.push_back("-distMap");
    opts.push_back("-go");
    printFormattedCommandLineOptions(commandLineOptions, "Transformation options", 100, opts);

    // Optimisation options
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-maxSwitchIt");
    opts.push_back("-ln");
    opts.push_back("-lp");
    opts.push_back("-maxFitIt");
    printFormattedCommandLineOptions(commandLineOptions, "Optimisation options", 100, opts);

    // Help/information
    opts.erase(opts.begin(), opts.end());
    opts.push_back("-h");
    printFormattedCommandLineOptions(commandLineOptions, "Help/information", 100, opts);

    return EXIT_SUCCESS;
  }

  //--------------------------------
  // Load the reference-state image 
  //--------------------------------
  /// \todo: invalid file name causes memory access violation
  std::string referenceStateImageFileName = parser->getCmdOptionAsString("-refState");
  nifti_image* referenceStateImage = nifti_image_read(referenceStateImageFileName.c_str(), true);

  if (referenceStateImage == nullptr)
  {
    char msg[200];
    sprintf_s(msg, "Could not read reference state image: %s", parser->getCmdOptionAsString("-refState").c_str());
    supremo_print_error(msg);
    supremo_exit(1, __FILE__, __LINE__);
  }

  //-------------------------
  // Load the dynamic images 
  //-------------------------

  int numberOfDynamicImages = parser->getCmdOptionAsInt("-dynamic", 0);
  std::string dynamicImageFileName = parser->getCmdOptionAsString("-dynamic", 1);
  nifti_image** dynamicImages = nullptr; // (nifti_image **)malloc( numberOfDynamicImages * sizeof( nifti_image * ) );
  std::ifstream dynamicNamesFile(dynamicImageFileName.c_str(), std::ifstream::in);

  if (!dynamicNamesFile.is_open())
  {
    char msg[200];
    sprintf_s(msg, "Cannot open the dynamic image names file %s", dynamicImageFileName.c_str());
    supremo_print_error(msg);
    supremo_exit(1, __FILE__, __LINE__);
  }

  // Read the complete file, check if correct number of files was specified in the file, then read
  {
    std::vector<std::string> allDynamicImageNames;
    std::string curDynamicImageName;

    // read values until no more a found in the file
    while (dynamicNamesFile >> curDynamicImageName)
    {
      allDynamicImageNames.push_back(curDynamicImageName);
    }

    // Check that the correct number of dynamic images was provided
    if (allDynamicImageNames.size() != numberOfDynamicImages)
    {
      supremo_print_error("Number of dynamic images not as expected.");
      supremo_exit(1, __FILE__, __LINE__);
    }

    // Now read the images
    dynamicImages = (nifti_image**)malloc(numberOfDynamicImages * sizeof(nifti_image*));

    for (int d = 0; d < numberOfDynamicImages; ++d)
    {
      dynamicImages[d] = nifti_image_read(allDynamicImageNames[d].c_str(), true);

      // Check if the current image was loaded properly
      if (dynamicImages[d] == nullptr)
      {
        char msg[200];
        sprintf_s(msg, "Unable to open dynamic image %i: %s", d, allDynamicImageNames[d].c_str());
        supremo_print_error(msg);
        supremo_exit(1, __FILE__, __LINE__);
      }
    }
  }
  dynamicNamesFile.close();


  //-------------------------
  // Load the surrogate data
  //-------------------------

  int numberOfSurrogateSignals = parser->getCmdOptionAsInt("-surr", 0);
  std::string surrogateFileName = parser->getCmdOptionAsString("-surr", 1);

  // Open the file  
  std::ifstream surrSignalFile(surrogateFileName.c_str(), std::ifstream::in);

  // Check file was opened correctly
  if (!surrSignalFile.is_open())
  {
    char msg[200];
    sprintf_s(msg, "Surrogate signal file could not be opened: %s", surrogateFileName.c_str());
    supremo_print_error(msg);
    supremo_exit(1, __FILE__, __LINE__);
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
      rawSurrValues.push_back(curSurrVal);
    }

    // Check that the correct number of surrogate signals was provided
    if (rawSurrValues.size() != numberOfSurrogateSignals * numberOfDynamicImages)
    {
      supremo_print_error("Number of surrogate signals not as expected.");
      supremo_exit(1, __FILE__, __LINE__);
    }

    // allocate float array of correct size and copy over values from the vector
    // we could just take a pointer like ( &rawSurrValues[0] ), but dynamic memory
    // management of the vector might invalidate that pointer. Hence explicit.
    surrogateSignals = new float[numberOfSurrogateSignals * numberOfDynamicImages];

    for (int i = 0; i < numberOfDynamicImages * numberOfSurrogateSignals; ++i)
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

  if (parser->cmdOptionExists("-defSpace"))
  {
    defSpaceImage = reg_io_ReadImageFile(parser->getCmdOptionAsString("-defSpace").c_str());

    if (defSpaceImage == nullptr)
    {
      char msg[200];
      sprintf_s(msg, "Could not read space definition image: %s", parser->getCmdOptionAsString("-defSpace").c_str());
      supremo_print_error(msg);
      supremo_exit(1, __FILE__, __LINE__);
    }
  }

  
  // Generate an image acquisition object as required
  std::shared_ptr<ImageAcquisition> imageAcquisition = nullptr;
  if (parser->cmdOptionExists("-dType"))
  {
    int iTmpDynType = parser->getCmdOptionAsInt("-dType");

    switch (iTmpDynType)
    {
    case SAME_RES_AS_STATIC:
      {
        auto noImageAcquisition = std::make_shared<NoImageAcquisition>();
        imageAcquisition = noImageAcquisition;
        break; 
      }
    case LOWER_RES_THAN_STATIC:
      {
        auto lowResImageAcquisition = std::make_shared<LowResolutionImageAcquisition>();
        imageAcquisition = lowResImageAcquisition;
        break;
      }
    default:
      {
        char msg[200];
        sprintf_s(msg, "Dynamic image type unknown: %i", iTmpDynType);
        supremo_print_error(msg);
        supremo_exit(1, __FILE__, __LINE__);
      }
    }
  }
  else
  {
    imageAcquisition = std::make_shared<NoImageAcquisition>();
  }



  // -mcrType
  // default to "no reconstruction" option
  // Build the reconstruction object here and pass it to supremo
  std::shared_ptr<MoCoRecon> motionCompensatedImageRecon = nullptr;

  // Motion compensated image reconstruction 
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
        avgWeightedReconstructor->SetImageAcquisition(imageAcquisition);
        motionCompensatedImageRecon = avgWeightedReconstructor;
        break;
      }
      case SUPER_RESOLUTION_RESTART:
      {
        auto superResReconstructor = std::make_shared<MoCoReconSuperResolutionIBP>( maxMCRIt, false );
        superResReconstructor->SetImageAcquisition(imageAcquisition);
        motionCompensatedImageRecon = superResReconstructor;
        break; 
      }
      case SUPER_RESOLUTION_UPDATE:
      {
        auto superResReconstructor = std::make_shared<MoCoReconSuperResolutionIBP>( maxMCRIt, true );
        superResReconstructor->SetImageAcquisition(imageAcquisition);
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
  supremo->SetImageAcquisition( imageAcquisition );
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
  std::size_t foundNiiInString = outputRCMFileName.find(".nii");
  if (foundNiiInString != std::string::npos)
  {
    outputRCMFileName.replace(foundNiiInString, outputRCMFileName.length(), "");
  }

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


  //----------
  // Clean up
  //----------

  // Delete the Supremo object
  supremo.reset();

  // Delete the command line parser
  parser.reset();

  // Delete the motion-compensated image reconstruction if applicable
  if (nullptr != motionCompensatedImageRecon)
  {
    motionCompensatedImageRecon.reset();
  }
  
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

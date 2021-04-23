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




#include "_reg_ReadWriteImage.h"
#include "Supremo.h"
#include "MoCoReconSuperResolutionIBP.h"
#include "SSDImageSimilarity.h"
#include "CommandLineParser.h"
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>



// Tolerances allowed 
constexpr auto EPS_SINGLE = 0.001;




int main( int argc, char **argv )
{
  // Note due to reduce computation time, here the image pyramid is actually calcualted (as opposed to the average weighted reconstruction)
  
  // Utilise the command line parser
  std::map<std::string, CommandLineOption> commandLineOptions;
  commandLineOptions["-refState"] = { 1, true, "Reference state image" };
  commandLineOptions["-dynamic"] = { 2, true, "Dynamic image data" };
  commandLineOptions["-surr"] = { 2, true, "Surrogate data" };
  commandLineOptions["-rcmIn"] = { 1, true, "Respiratory correspondence model in" };
  commandLineOptions["-mcrToCompare"] = { 1, true, "The reference motion-compensated reconstruction." };
  commandLineOptions["-ln"] = { 1, true, "Number of levels performed in the original MCR." };
  commandLineOptions["-mcrIt"] = { 1, true, "Number of iterations performed during reconstruction." };
  commandLineOptions["-outMCR"] = { 1, false, "Motion-compensated image - for visual comparison." };

  // Parse the command line
  std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>( argc, argv, commandLineOptions );
  std::cout << parser->getCommandLine() << std::endl;

  const unsigned int numberOfLevelsToPerform = static_cast<unsigned int>(parser->getCmdOptionAsInt( "-ln" ));
  const unsigned int mcrIterations = static_cast<unsigned int>(parser->getCmdOptionAsInt( "-mcrIt" ));

  // Read the reference state image
  const std::string referenceStateImageFileName = parser->getCmdOptionAsString( "-refState" );
  nifti_image* referenceStateImage = reg_io_ReadImageFile( referenceStateImageFileName.c_str() );

  // Read the correspondence model image
  std::string rcmImageFileName = parser->getCmdOptionAsString( "-rcmIn" );
  nifti_image* rcmImage = reg_io_ReadImageFile( rcmImageFileName.c_str() );

  // Read the expected reconstructed image
  std::string mcrImageFileName = parser->getCmdOptionAsString( "-mcrToCompare" );
  nifti_image* mcrImage = reg_io_ReadImageFile( mcrImageFileName.c_str() );


  if (referenceStateImage == nullptr)
  {
    char msg[200];
    sprintf_s( msg, "Could not read reference state image: %s", parser->getCmdOptionAsString( "-refState" ).c_str() );
    supremo_print_error( msg );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  //-------------------------
  // Load the dynamic images 
  //-------------------------
  int numberOfDynamicImages = parser->getCmdOptionAsInt( "-dynamic", 0 );
  std::string dynamicImageFileName = parser->getCmdOptionAsString( "-dynamic", 1 );
  std::vector<nifti_image*> dynamicImages;
  std::ifstream dynamicNamesFile( dynamicImageFileName.c_str(), std::ifstream::in );

  if (!dynamicNamesFile.is_open())
  {
    char msg[200];
    sprintf_s( msg, "Cannot open the dynamic image names file %s", dynamicImageFileName.c_str() );
    supremo_print_error( msg );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
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
      supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
    }

    // Now read the images
    for (int d = 0; d < numberOfDynamicImages; ++d)
    {
      dynamicImages.push_back( reg_io_ReadImageFile( allDynamicImageNames[d].c_str() ) );

      // Check if the current image was loaded properly
      if (dynamicImages[d] == nullptr)
      {
        char msg[200];
        sprintf_s( msg, "Unable to open dynamic image %i: %s", d, allDynamicImageNames[d].c_str() );
        supremo_print_error( msg );
        supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
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
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
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
      supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
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

  std::vector<std::vector<float>> surrogateSignalsVec;
  // convert to vector
  for (int iTimePoint = 0; iTimePoint < numberOfDynamicImages; ++iTimePoint)
  {
    std::vector<float> tmpVect;
    for (int iSurrSig = 0; iSurrSig < numberOfSurrogateSignals; ++iSurrSig)
    {
      tmpVect.push_back( surrogateSignals[iTimePoint * numberOfSurrogateSignals + iSurrSig] );
    }
    surrogateSignalsVec.push_back( tmpVect );
  }

  // Generate the image pyramid
  std::vector<std::shared_ptr<ImagePyramid<float> > > dynamicPyramids;
  std::shared_ptr<ImagePyramid<float> > referenceStatePyramid;

  // Generate the dynamic image pyramids
  for (unsigned int i = 0; i < numberOfDynamicImages; ++i)
  {
    auto curDynPyramid = std::make_shared<ImagePyramid<float> >();
    curDynPyramid->GenerateLevels( dynamicImages[i], numberOfLevelsToPerform, numberOfLevelsToPerform );
    dynamicPyramids.push_back( curDynPyramid );
  }
  
  // Generate the reference state image pyramid
  referenceStatePyramid = std::make_shared<ImagePyramid<float> >();
  referenceStatePyramid->GenerateLevels( referenceStateImage, numberOfLevelsToPerform, numberOfLevelsToPerform );

  // Extract relevant dynamic images from  vector of pyramids
  std::vector<nifti_image*> curDynamicImages;
  for (int nDynImg = 0; nDynImg < numberOfDynamicImages; ++nDynImg)
  {
    curDynamicImages.push_back( dynamicPyramids[nDynImg]->GetLevel( 0 ) );
  }

  // Set up related classes
  // Transformation
  float bSplineCPGSpacing[3] = { rcmImage->dx, rcmImage->dy, rcmImage->dz };
  auto bsplTrafo = std::make_shared<BSplineTransformation>( referenceStateImage, 1, bSplineCPGSpacing );
  bsplTrafo->InitialiseLevel( 0 );

  // Correspondnece model
  auto correspModel = std::make_shared<CorrespondenceModel>( numberOfSurrogateSignals, bsplTrafo );

  // Only initialise the first level
  correspModel->InitialiseLevel( 0 );
  
  // Feed the input data into the correspondence model
  // but before check that the data is of expected type
  if (rcmImage->nbyper != sizeof( CorrespondenceModel::PrecisionType ))
  {
    supremo_print_error( "Data-type mismatch between input RCM and correspondence model class." );
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }
  correspModel->SetParameters( (CorrespondenceModel::PrecisionType*) rcmImage->data );

  // image acquisition
  std::shared_ptr<ImageAcquisition> lowResAcquisition = std::make_shared<LowResolutionImageAcquisition>();

  // Motion-compensated image reconstruction to test
  auto moCoReco = std::make_shared<MoCoReconSuperResolutionIBP>( mcrIterations, false );
  moCoReco->SetCorrespondenceModel( correspModel );
  moCoReco->SetReconstructionGeometryImage( referenceStatePyramid->GetLevel(0) );
  moCoReco->SetSurrogateSignals( surrogateSignalsVec );
  moCoReco->SetDynamicImages( curDynamicImages );
  moCoReco->SetImageAcquisition( lowResAcquisition );
  moCoReco->Update();

  // Save the reconstructed image if requried
  if (parser->cmdOptionExists( "-outMCR" ))
  {
    reg_io_WriteImageFile( moCoReco->GetReconstructedImage(), parser->getCmdOptionAsString( "-outMCR" ).c_str() );
  }

  // Compare the output with the expected outcome
  return checkImageEquality<float>( moCoReco->GetReconstructedImage(), mcrImage, EPS_SINGLE, 2 );
}

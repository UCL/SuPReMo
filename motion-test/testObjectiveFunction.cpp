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
#include "SSDImageSimilarity.h"
#include "CommandLineParser.h"
#include "ObjectiveFunction.h"
#include "NoImageAcquisition.h"
#include <memory>
#include <cmath>
#include <fstream>



// Tolerances allowed 
constexpr auto EPS_SINGLE = 0.0001;
constexpr auto EPS_SINGLE_ABS = 0.001;
constexpr auto EPS_SINGLE_REL = 0.001;

bool AlmostEqualRelativeAndAbs(float A, float B,
	float maxDiff, float maxRelDiff = FLT_EPSILON)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	float diff = fabs(A - B);
	if (diff <= maxDiff)
		return true;

	A = fabs(A);
	B = fabs(B);
	float largest = (B > A) ? B : A;

	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}




int main(int argc, char **argv)
{

	// Some fixed values for now
	double bSplineBendingEnergyWeight = 0.1;
	double bSplineLinearEnergyWeight = 0.0;
	unsigned int numberOfLevelsToPerform = 1;
	float bSplineCPGSpacing[3] = { 10.,10.,10. };

	// Utilise the command line parser
	std::map<std::string, CommandLineOption> commandLineOptions;
	commandLineOptions["-refState"] = { 1, true, "Reference state image" };
	commandLineOptions["-dynamic"] = { 2, true, "Dynamic image data" };
	commandLineOptions["-surr"] = { 2, true, "Surrogate data" };
	commandLineOptions["-valFile"] = { 1, true, "Expected outcome at zero" };
	commandLineOptions["-gradFiles"] = { 1, true, "Expected gradient at zero (as nifti image), komma separated list." };
	commandLineOptions["-oGradImgs"] = { 1, false, "File names of the output gradient images. " };
	// Parse the command line
	std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>(argc, argv, commandLineOptions);
	std::cout << parser->getCommandLine() << std::endl;

	std::string referenceStateImageFileName = parser->getCmdOptionAsString("-refState");
	nifti_image* referenceStateImage = reg_io_ReadImageFile(referenceStateImageFileName.c_str());

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
		dynamicImages = (nifti_image **)malloc(numberOfDynamicImages * sizeof(nifti_image *));

		for (int d = 0; d < numberOfDynamicImages; ++d)
		{
			dynamicImages[d] = reg_io_ReadImageFile(allDynamicImageNames[d].c_str());

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
			tmpVect.push_back(surrogateSignals[iTimePoint * numberOfSurrogateSignals + iSurrSig]);
		}
		surrogateSignalsVec.push_back(tmpVect);
	}

	// Load the gradient images
	auto gradFileNames = splitStringbyDelimiter(parser->getCmdOptionAsString("-gradFiles"), ",");
	if (gradFileNames.size() != numberOfSurrogateSignals)
	{
		supremo_print_error("Number of gradient images not as expected.");
		supremo_exit(1, __FILE__, __LINE__);
	}

	// Read the gradient images and get the number of parameters
	size_t numOfParametersInputGrad = 0;
	std::vector<nifti_image*> gradImages;
	for (int iSurrSignals = 0; iSurrSignals < numberOfSurrogateSignals; ++iSurrSignals)
	{
		nifti_image* curGradImg = reg_io_ReadImageFile(gradFileNames[iSurrSignals].c_str());
		if (curGradImg == NULL)
		{
			supremo_print_error("The gradient image could not be read");
			return EXIT_FAILURE;
		}
		gradImages.push_back(curGradImg);
		numOfParametersInputGrad += curGradImg->nvox;
	}

	// Set up related classes
	auto bsplTrafo = std::make_shared<BSplineTransformation>(referenceStateImage,
		numberOfLevelsToPerform,
		bSplineCPGSpacing);
	bsplTrafo->SetBendingEnergyWeight(bSplineBendingEnergyWeight);
	bsplTrafo->SetLinearEnergyWeight(bSplineLinearEnergyWeight);
	bsplTrafo->InitialiseLevel(0);

	std::shared_ptr<SSDImageSimilarity> simMeasure = std::make_shared<SSDImageSimilarity>();
	std::shared_ptr<CorrespondenceModel> correspondenceModel = std::make_shared<CorrespondenceModel>(numberOfSurrogateSignals, bsplTrafo);
	correspondenceModel->InitialiseLevel(0);

  std::shared_ptr<ImageAcquisition> imageAcquisition = std::make_shared<NoImageAcquisition>();

	std::shared_ptr<ObjectiveFunction> objectiveFunction = std::make_shared<ObjectiveFunction>();
  objectiveFunction->SetCorrespondenceModel( correspondenceModel );
  objectiveFunction->SetSurrogateSignals( surrogateSignalsVec );
  objectiveFunction->SetSimilarityMeasure( simMeasure );
  objectiveFunction->SetReferenceStateImage( referenceStateImage );
  objectiveFunction->SetImageAcquisition( imageAcquisition );
	std::vector<nifti_image*> curDynamicImages;
	for (int nDynImg = 0; nDynImg < numberOfDynamicImages; ++nDynImg)
	{
    curDynamicImages.push_back( dynamicImages[nDynImg] );
	}
	objectiveFunction->SetDynamicImages(curDynamicImages);

  CorrespondenceModel::PrecisionType* corParams = (CorrespondenceModel::PrecisionType*)calloc( correspondenceModel->GetNumberOfParameters(), sizeof( CorrespondenceModel::PrecisionType ) );
  CorrespondenceModel::PrecisionType* gradAtZero = (CorrespondenceModel::PrecisionType*)calloc( correspondenceModel->GetNumberOfParameters(), sizeof( CorrespondenceModel::PrecisionType ) );
  CorrespondenceModel::PrecisionType measuredValueAtZero = objectiveFunction->GetValue( corParams );
  objectiveFunction->GetGradient( corParams, gradAtZero );


	// Open the expected outcome file
	std::ifstream objFuncValAtZeroFile(parser->getCmdOptionAsString("-valFile").c_str(), std::ifstream::in);

	if (!objFuncValAtZeroFile.is_open())
	{
		char msg[200];
		sprintf_s(msg, "Fils with the expected value could not be obened %s", surrogateFileName.c_str());
		supremo_print_error(msg);
    supremo_exit( 1, __FILE__, __LINE__ );
	}
	float expectedObjFuncValAtZero;
	objFuncValAtZeroFile >> expectedObjFuncValAtZero;

  // Make sure that all tests complete
  bool testSuccessful = true;
  
  if (fabs(measuredValueAtZero - expectedObjFuncValAtZero) > EPS_SINGLE)
	{
    char msg[400];
    sprintf_s(msg, "Measured objective function value at zero not as expected. Evaluated: %0.6f, expected: %0.6f", measuredValueAtZero, expectedObjFuncValAtZero );
    supremo_print_error( msg );
    testSuccessful = false;
	}
  
	// Chech the number of parameters
	if (numOfParametersInputGrad != objectiveFunction->GetNumberOfParameters())
	{
		supremo_print_error("Number of objective function parameters not as expected.");
    testSuccessful = false;
	}

	// Generate output gradient images to save as output if requried
	if (parser->cmdOptionExists("-oGradImgs"))
	{
		auto outFileNames = splitStringbyDelimiter(parser->getCmdOptionAsString("-oGradImgs"), ",");
		
		for (size_t iImgCount = 0; iImgCount < gradImages.size(); ++iImgCount)
		{
			size_t offset = iImgCount * gradImages[iImgCount]->nvox;
			nifti_image* curImg = nifti_copy_nim_info(gradImages[iImgCount]);
			curImg->data = malloc(gradImages[iImgCount]->nvox * gradImages[iImgCount]->nbyper);
			memcpy(curImg->data, &(gradAtZero[offset]), gradImages[iImgCount]->nvox * gradImages[iImgCount]->nbyper);
			reg_io_WriteImageFile(curImg, outFileNames[iImgCount].c_str());
			free(curImg->data);
			curImg->data = nullptr;
			nifti_image_free(curImg);
		}
	}

	// Check the gradient values
	unsigned int iContCount = 0;
	for (size_t iImgCount = 0; iImgCount < gradImages.size(); ++iImgCount )
	{
		float* pCurImgData = (float*) gradImages[iImgCount]->data;
		
		for (size_t iVoxCount = 0; iVoxCount < gradImages[iImgCount]->nvox; ++iVoxCount)
		{
			if (!AlmostEqualRelativeAndAbs(gradAtZero[iContCount], pCurImgData[iVoxCount], EPS_SINGLE_ABS, EPS_SINGLE_REL ) )
			{
				supremo_print_error("Gradient value not as expected.");
        testSuccessful = false;
			}
			iContCount++;
		}	
	}


  if (testSuccessful)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

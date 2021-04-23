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
#include "CommandLineParser.h"
#include "checkEquality.h"
#include <memory>
#include <cmath>
#include <fstream>


// Tolerances allowed 
constexpr auto EPS_SINGLE = 0.0001;
constexpr auto EPS_SINGLE_ABS = 0.001;
constexpr auto EPS_SINGLE_REL = 0.001;




int main(int argc, char **argv)
{
	// Utilise the command line parser
	std::map<std::string, CommandLineOption> commandLineOptions;
	commandLineOptions["-cpg"]    = { 1, true, "Control-point grid image" };
	commandLineOptions["-cpgRef"] = { 1, true, "Reference image that was used to generate the above CPG." };
	commandLineOptions["-imgX"]   = { 1, true, "First image to be transformed" };
	commandLineOptions["-imgY"]   = { 1, true, "Second image to be transformed with the adjoint" };
	
  // Parse the command line
	std::shared_ptr<CommandLineParser> parser = std::make_shared<CommandLineParser>(argc, argv, commandLineOptions);
	std::cout << parser->getCommandLine() << std::endl;

  // Read the CPG image
  std::string cpgImageFileName = parser->getCmdOptionAsString( "-cpg" );
  nifti_image* cpgImage = reg_io_ReadImageFile( cpgImageFileName.c_str() );
  if (cpgImage == NULL)
  {
    supremo_print_error( "The cpg image could not be read" );
    return EXIT_FAILURE;
  }

  // Read the reference image used for the CPG generation
  std::string cpgReferenceImageFileName = parser->getCmdOptionAsString( "-cpgRef" );
  nifti_image *cpgReferenceImage = reg_io_ReadImageFile( cpgReferenceImageFileName.c_str() );
  if (cpgReferenceImage == NULL)
  {
    supremo_print_error( "The cpg reference image could not be read" );
    return EXIT_FAILURE;
  }

  // Read imgX
  std::string imageXFileName = parser->getCmdOptionAsString( "-imgX" );
  nifti_image* imageX = reg_io_ReadImageFile( imageXFileName.c_str() );
  if (imageX == NULL || imageX->datatype != DT_FLOAT32)
  {
    supremo_print_error( "The image X could not be read or is not of type float" );
    return EXIT_FAILURE;
  }

  // Read imgY  
  std::string imageYFileName = parser->getCmdOptionAsString( "-imgY" );
  nifti_image* imageY = reg_io_ReadImageFile( imageYFileName.c_str() );
  if (imageY == NULL || imageY->datatype != DT_FLOAT32)
  {
    supremo_print_error( "The image Y could not be read or is not of type float" );
    return EXIT_FAILURE;
  }

  // Get the cpg spacing from the image information
  float* cpgSpacing = new float[cpgImage->ndim];

  if (cpgImage->ndim >= 2)
  {
    cpgSpacing[0] = cpgImage->dx;
    cpgSpacing[1] = cpgImage->dy;
  }

  if (cpgImage->ndim > 2)
  {
    cpgSpacing[2] = cpgImage->dz;
  }

  // Set up the transformation
  // Generate the b-spline transformation object
  std::shared_ptr<BSplineTransformation> bsplTrafo = std::make_shared<BSplineTransformation>( cpgReferenceImage, 1, cpgSpacing );

  // Check that the generated B-spline transformation and the input cpg have the same number of parameters
  if (bsplTrafo->GetNumberOfParameters() != cpgImage->nvox)
  {
    supremo_print_error( "Wrong number of parameters in BSpline-Transformation" );
    return EXIT_FAILURE;
  }

  // Force the CPG values into the transformation object
  bsplTrafo->SetParameters( (float*)cpgImage->data, false );
  
  // Calculate Tx has to have the shape of y, thus select this as the reference image
  nifti_image* imageXTransformed = nifti_copy_nim_info( imageY );
  imageXTransformed->data = (void*)calloc( imageXTransformed->nvox, imageXTransformed->nbyper );
  bsplTrafo->TransformImage( imageX, imageXTransformed );


  // Calculate T* y
  
  // Prepare the output and the weights
  nifti_image* imageYWeights             = nifti_copy_nim_info( imageY );
  nifti_image* imageYAfterAdjoint        = nifti_copy_nim_info( imageY );
  nifti_image* imageYAfterAdjointWeights = nifti_copy_nim_info( imageY );

  imageYWeights->data             = calloc( imageYWeights->nvox,             imageYWeights->nbyper             );
  imageYAfterAdjoint->data        = calloc( imageYAfterAdjoint->nvox,        imageYAfterAdjoint->nbyper        );
  imageYAfterAdjointWeights->data = calloc( imageYAfterAdjointWeights->nvox, imageYAfterAdjointWeights->nbyper );

  // fill the source weights with appropriate values
  float* imageYWeightsPointer = static_cast<float*>(imageYWeights->data);
  for (unsigned int i = 0; i < imageYWeights->nvox; ++i)
  {
    imageYWeightsPointer[i] = 1.f;
  }


  // Note: the DVF is calcualted for the size of the source image internally
  bsplTrafo->TransformImageAdjoint( imageY, imageYWeights, imageYAfterAdjoint, imageYAfterAdjointWeights );

  // Normalise after transformation with adjoint
  float* adjImgPtr = static_cast<float*>(imageYAfterAdjoint->data);
  float* adjWeigthsPtr = static_cast<float*>(imageYAfterAdjointWeights->data);
  
  // Note: This normalisation is required to NOT make intensity changes according to 
  //       diverging or converging DVFs. Such intensity changes however are an essential 
  //       assumption of the tests based on the inner product. To fully comply with this 
  //       assumption the pull interpolation should also implement intensity changes based
  //       on the Jacobian determinant - which is currently not the case and will need 
  //       to be looked at in the future. 
  //for (size_t i = 0; i < imageYAfterAdjoint->nvox; ++i)
  //{
  //  if (adjWeigthsPtr[i] != 0)
  //    adjImgPtr[i] = adjImgPtr[i] / adjWeigthsPtr[i];
  //  else
  //    adjImgPtr[i] = 0.f;
  //}
  

  // Calculate the inner products
  // < Tx ,y > and < x, T* y > 
  
  size_t numVox = imageYAfterAdjoint->nvox;
  if (numVox != imageXTransformed->nvox || numVox != imageX->nvox || numVox != imageY->nvox)
  {
    supremo_print_error( "The image sizes generated for comparison do not match..." );
    return EXIT_FAILURE;
  }
  
  float innerProdTxy  = 0.f;
  float innerProdxTsy = 0.f;
  float* ptrX   = static_cast<float*>(imageX->data);
  float* ptrTx  = static_cast<float*>(imageXTransformed->data);
  float* ptrY   = static_cast<float*>(imageY->data);
  float* ptrTsy = static_cast<float*>(imageYAfterAdjoint->data);
  
  for (size_t i = 0; i < numVox; ++i)
  {
    if ((!isnan( ptrTx[i] )) && (!isnan( ptrY[i] )))
    {
      innerProdTxy += ptrTx[i] * ptrY[i];
    }
    
    if ((!isnan( ptrX[i] )) && (!isnan( ptrTsy[i] )))
    {
      innerProdxTsy += ptrX[i] * ptrTsy[i];
    }
  }

  std::cout << " <Tx , y    > = " << innerProdTxy << std::endl;
  std::cout << " < x , T* y > = " << innerProdxTsy << std::endl;
  
  // comapre
  // Relative difference, use scaling according to image size
  innerProdTxy  /= (float)numVox;
  innerProdxTsy /= (float)numVox;
  
  float relDiff = abs(innerProdTxy - innerProdxTsy) / abs(innerProdxTsy);
  std::cout << "|<Tx,y> - <x,T*y>| / |<x,T*y>| = " << relDiff << std::endl;

  if ( relDiff > EPS_SINGLE_ABS)
  {
    std::cout << "Dfference between inner products exceeds allowed limits." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

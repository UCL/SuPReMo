#include "CorrespondenceModel.h"
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
#include "CorrespondenceModel.h"




//------------------------------------------
// CorrespondenceModel::CorrespondenceModel
//------------------------------------------
CorrespondenceModel::CorrespondenceModel( unsigned int numberOfSurrogateSignalsIn,
  std::shared_ptr<Transformation> transformationIn ) :
  numberOfSurrogateSignals( numberOfSurrogateSignalsIn ),
  numberOfTransformationParameters( 0 ),
  transform( nullptr ),
  modelParameters( nullptr ),
  recoveryModelParameters( nullptr ),
  numberOfModelParameters( 0 )
{
  this->transform = transformationIn->DeepCopy();
  this->numberOfTransformationParameters = this->transform->GetNumberOfParameters();
  this->numberOfModelParameters = this->numberOfTransformationParameters * this->numberOfSurrogateSignals;

  // Allocate internal model parameter array
  // Should be contiguous thus not a separate one per surrogate signal
  this->modelParameters = new PrecisionType[this->numberOfModelParameters]();

  if (this->modelParameters == nullptr)
  {
    char msg[200];
    sprintf_s( msg, "Memory allocation for correspondence model parameters failed." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
}




//-------------------------------------------
// CorrespondenceModel::~CorrespondenceModel
//-------------------------------------------
CorrespondenceModel::~CorrespondenceModel()
{
  // Remove the reference to the transformation pointer
  this->transform = nullptr;
  
  // Delete the model parameters
  if ( this->modelParameters != nullptr )
  {
    delete[] this->modelParameters;
    this->modelParameters = nullptr;
  }

  // Delete the recovery parameters if they were used
  this-> ClearRecoveryModelParameters();
}




//-----------------------------------------------------------
// CorrespondenceModel::GetTransformationFromSurrogateSignal
//-----------------------------------------------------------
std::shared_ptr<Transformation> CorrespondenceModel::GetTransformationFromSurrogateSignal( const SurrogateSignalType & surrogateSignalIn )
{
  // Overview
  // 0) Check the size of the input vector
  // 1) Set the parameters of the internal transformation to the requested value
  // 2) Generate a (deep) copy of the transformation and return it 
  //    ToDo: Need to review if a copy is actually required. 
  if ( surrogateSignalIn.size() != this->numberOfSurrogateSignals )
  {
    // Exit    
    char msg[200];
    sprintf_s( msg, "Surrogate signal of wrong dimension provided." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Set the parameters further down, to enable conversion from displacement convention to deformation.
  Transformation::PrecisionType* pTransformParameters = (Transformation::PrecisionType*)malloc( this->transform->GetNumberOfParameters()*sizeof( Transformation::PrecisionType ) );

  if ( nullptr == pTransformParameters )
  { 
    supremo_print_error("Could not allocate memory");
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // indices to access the surrogate signal and the transformation parameters
  unsigned int iSurSig = 0;
  PrecisionType surrSigValue = surrogateSignalIn[iSurSig];

  for ( unsigned int iTransformParam = 0; iTransformParam < this->numberOfTransformationParameters; ++iTransformParam )
  {
    // Set the values in the first round
    pTransformParameters[iTransformParam] = surrSigValue * this->modelParameters[iTransformParam];
  }
  for ( iSurSig=1; iSurSig < this->numberOfSurrogateSignals; ++iSurSig )
  {
    surrSigValue = surrogateSignalIn[iSurSig];

    for ( unsigned int iTransformParam = 0; iTransformParam < this->numberOfTransformationParameters; ++iTransformParam )
    {
      // Add to the parameters in the second and further rounds
      pTransformParameters[iTransformParam] += surrSigValue * this->modelParameters[iTransformParam + iSurSig*this->numberOfTransformationParameters];
    }
  }
  
  this->transform->SetParameters( pTransformParameters, true ); 

  free( pTransformParameters );

  //return this->transform->DeepCopy();
  return this->transform;
}




//--------------------------------------
// CorrespondenceModel::InitialiseLevel
//--------------------------------------
void CorrespondenceModel::InitialiseLevel( unsigned int levelIn )
{
  // Overview:
  // 1) For each surrogate signal set the corresponding model parameters to the internal transform
  // 2) Generate a deep copy of this transform
  // 3) Refine the transformation parameters
  // 4) Allocate memory for all model parameters and copy the refined transformation parameters
  // 5) Delete the recovery parameters 

  // Vector holding temporary copies of the refined transformations. One for each surrogate signal
  std::vector<std::shared_ptr<Transformation> > refinedTransforms;

  // Allocate space for the transformation parameters
  //Transformation::PrecisionType* pTransformParameters = this->transform->GetParameters();
  Transformation::PrecisionType* pTransformParameters  = new Transformation::PrecisionType[this->transform->GetNumberOfParameters()]();

  for ( unsigned int iSurSig = 0; iSurSig < this->numberOfSurrogateSignals; ++iSurSig )
  {
    // Copy the model parameters corresponding to the current surrogate signal into the transformation
    for ( unsigned int iTrafoParams = 0; iTrafoParams < this->numberOfTransformationParameters; ++iTrafoParams )
    {
      pTransformParameters[iTrafoParams] = this->modelParameters[iTrafoParams + iSurSig * this->numberOfTransformationParameters];
    }
    
    // Note: The parameters actually are displacements despite setting the sceond parameter to false. However, 
    //       we just use the transformation to upscale the displacements, so on setting the parameters, we do not want to 
    //       convert to deformations
    this->transform->SetParameters( pTransformParameters, false );
    
    // Generate a copy with these parameters 
    refinedTransforms.push_back( this->transform->DeepCopy() );
  }
  // Delete the memory allocated for generating the intermediate transforms
  delete[] pTransformParameters;

  // Refine all copied transformations
  for ( unsigned int iSurSig = 0; iSurSig < this->numberOfSurrogateSignals; ++iSurSig )
  {
    refinedTransforms[iSurSig]->InitialiseLevel( levelIn );
  }

  // Update the number of parameters of the transformation and the correspondence model
  this->numberOfTransformationParameters = refinedTransforms[0]->GetNumberOfParameters();
  this->numberOfModelParameters = this->numberOfTransformationParameters * this->numberOfSurrogateSignals;

  // Free the old memory and allocate one for the new size
  delete[] this->modelParameters;
  this->modelParameters = new PrecisionType[ this->numberOfModelParameters ]();

  // Copy the upscaled model parameters
  for ( unsigned int iSurSignal = 0; iSurSignal < this->numberOfSurrogateSignals; ++iSurSignal )
  {
    auto curTrafoParams = refinedTransforms[iSurSignal]->GetCopyOfParameters();
    
    for ( unsigned int iTrafoParam = 0; iTrafoParam < this->numberOfTransformationParameters; ++iTrafoParam )
    {
      this->modelParameters[iTrafoParam + iSurSignal*this->numberOfTransformationParameters] = curTrafoParams[iTrafoParam];
    }
    free( curTrafoParams );
  }

  // Update the internal transformation
  this->transform.reset();
  this->transform = refinedTransforms[0];

  // Delete the recovery parameters if they were used
  // memory must not be used on a different resolution level
  if (nullptr != this->recoveryModelParameters)
  {
    delete[] this->recoveryModelParameters;
    this->recoveryModelParameters = nullptr;
  }
}




//----------------------------------------
// CorrespondenceModel::GetTransformation
//----------------------------------------
std::shared_ptr<Transformation>  CorrespondenceModel::GetTransformation()
{
  // Simply return the shared pointer to the transformation of this correspondence model
  return this->transform;
}




//----------------------------------------------------
// CorrespondenceModel::GetCorrespondenceModelAsImage
//----------------------------------------------------
std::vector<nifti_image*> CorrespondenceModel::GetCorrespondenceModelAsImage()
{
	// Get a structure ready to save the transformations
	std::vector<std::shared_ptr<Transformation> > vTransforms;
	
	// Generate a transformation object for every surrogate signal
	for (unsigned int iSurrSig = 0; iSurrSig < this->numberOfSurrogateSignals; ++iSurrSig)
	{
		PrecisionType* ptrFirstParam = &(this->modelParameters[iSurrSig * this->numberOfTransformationParameters]);
		
		// Make sure that the displacements are not converted to deformations, hence use parameters as they are (i.e. called with false)
		this->transform->SetParameters( ptrFirstParam, false );
    vTransforms.push_back( this->transform->DeepCopy() );
	}

  // Each transformation may return multiple images (e.g. for sliding regions). These images do not have to be of the same size or geometry etc.
  // But each surrogate signal has the same transformation (deep copy above) and thus number of parameters

	// Generate the output images
  std::vector<nifti_image*> returnedCorrespondenceModelImages;

  // Here we need to know the number of images per transformation
  std::vector<nifti_image*> firstTrafoImgs = vTransforms[0]->GetTransformationAsImage();

  // Iterate over all images of the transformation and generate a correspondence model image for each transformation image
  // The correspondence model images will have a higher dimensionality, i.e. dim[6] will have the size of the number of surrogate signals
  for (size_t nTrafoImg = 0; nTrafoImg < firstTrafoImgs.size(); ++nTrafoImg)
  {
    // Allocate the correspondence model image header (only need to access the first transformation since the others are the same) ...
    nifti_image* curCorrespondenceModelImage = nifti_copy_nim_info( firstTrafoImgs[nTrafoImg] );
    
    // Modify the nifti header to accommodate the higher dimensionality of the correspondence model
    curCorrespondenceModelImage->dim[0] = curCorrespondenceModelImage->ndim = 6;
    curCorrespondenceModelImage->dim[6] = curCorrespondenceModelImage->nv = this->numberOfSurrogateSignals;
    curCorrespondenceModelImage->nvox = curCorrespondenceModelImage->nx *
      curCorrespondenceModelImage->ny *
      curCorrespondenceModelImage->nz *
      curCorrespondenceModelImage->nt *
      curCorrespondenceModelImage->nu *
      curCorrespondenceModelImage->nv;
    
    // Allocate the data and get a pointer to this
    curCorrespondenceModelImage->data = (void*)malloc( curCorrespondenceModelImage->nvox * curCorrespondenceModelImage->nbyper );
    PrecisionType *curCorrespModDataPtr = static_cast<PrecisionType *>(curCorrespondenceModelImage->data);

    // Now iterate over all surrogate signals and fill the previously allocate image
    for (unsigned int nSurSig = 0; nSurSig < this->numberOfSurrogateSignals; nSurSig++)
    {
      auto curTransformImg = vTransforms[nSurSig]->GetTransformationAsImage()[nTrafoImg];

      memcpy( &curCorrespModDataPtr[nSurSig * curTransformImg->nvox],
        curTransformImg->data, curTransformImg->nvox * curTransformImg->nbyper );
    }
    
    // This higher dimensional image will be returned.
    returnedCorrespondenceModelImages.push_back( curCorrespondenceModelImage );
  }

  return returnedCorrespondenceModelImages;
}




//------------------------------------------------------------
// CorrespondenceModel::SaveCurrentModelParametersForRecovery
//------------------------------------------------------------
void CorrespondenceModel::SaveCurrentModelParametersForRecovery()
{
  // Allocate the memory if it does not already exist
  // This memory will be deleted if a new resolution level is initialised, so it does not survive
  // the resolution level switch. 
  if (nullptr == this->recoveryModelParameters) 
  {
    this->recoveryModelParameters = new PrecisionType[this->numberOfModelParameters]();
  }
  
  // Copy over the memory
  memcpy( this->recoveryModelParameters, this->modelParameters, sizeof( PrecisionType )*this->numberOfModelParameters );
}




//--------------------------------------------------
// CorrespondenceModel::RecoverSavedModelParameters
//--------------------------------------------------
void CorrespondenceModel::RecoverSavedModelParameters()
{
  // Make sure that model parameters were previously saved for recovery
  if (nullptr == this->recoveryModelParameters)
  {
    supremo_print_error("No model parameters saved to recover");
    supremo_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Copy the stored model parameters to the current one
  memcpy( this->modelParameters, this->recoveryModelParameters, sizeof( PrecisionType ) * this->numberOfModelParameters );
  
  // And clean up the recovery memory
  this->ClearRecoveryModelParameters();
}




//---------------------------------------------------
// CorrespondenceModel::ClearRecoveryModelParameters
//---------------------------------------------------
void CorrespondenceModel::ClearRecoveryModelParameters()
{
  // Delete the recovery parameters
  if (nullptr != this->recoveryModelParameters)
  {
    delete[] this->recoveryModelParameters;
    this->recoveryModelParameters = nullptr;
  }
}




//------------------------------------
// CorrespondenceModel::GetParameters
//------------------------------------
CorrespondenceModel::PrecisionType * CorrespondenceModel::GetParameters()
{
  // Return the pointer to the model parameters
  return this->modelParameters;
}




//------------------------------------
// CorrespondenceModel::SetParameters
//------------------------------------
void  CorrespondenceModel::SetParameters( const CorrespondenceModel::PrecisionType * parametersIn)
{
  // Copy the parameters to the internal model parameters
  memcpy( this->modelParameters, parametersIn, sizeof( PrecisionType ) * this->numberOfModelParameters );
}




//------------------------------------
// CorrespondenceModel::SetParameters
//------------------------------------
void CorrespondenceModel::SetParameters( const std::vector<nifti_image*>& parameterImagesIn )
{
  // The images are arranged such that we have a correspondence model image "per sliding region"
  // The transformation internally expects the data as follows, where RX are the parameters per region:
  //  P_T = | R1 | R2 |  .....  | RN |
  //
  // This correspondence model concatenates the transformation parameters for each surrogate signal (M surrogate signals). 
  // Each transformation has by definition the same number of parameters 
  //  T1 = | R1 | R2 |  .....  | RN |
  //  T2 = | R1 | R2 |  .....  | RN |
  //  ...
  //  TM = | R1 | R2 |  .....  | RN |
  
  // The images given as an input are for each RX with dimension nv = M.
  //       | R1 | R2 |  .....  | RN |
  //       | R1 | R2 |  .....  | RN |
  //  
  //       | R1 | R2 |  .....  | RN |
  //        ==== ====           ====
  //        Img1 Img2           Img3
  
  // Check that the input fits to the current parameter configuration as described above

  // Check the precision type first
  for (size_t nImg = 0; nImg < parameterImagesIn.size(); ++nImg)
  {
    if (parameterImagesIn[nImg]->nbyper != sizeof( PrecisionType ))
    {
      supremo_print_error( "Precision type does not match. " );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }
  
  // The number of surrogate signals must be the same as the dimension nv in the given images
  for (size_t nImg = 0; nImg < parameterImagesIn.size(); ++nImg)
  {
    if (parameterImagesIn[nImg]->nv != this->numberOfSurrogateSignals)
    {
      supremo_print_error( "Number of surrogate signals does not correspond to nv dimension of input image." );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
  }

  // The number of voxels for every "nv" dimension must be the number of transformation parameters
  // Thus testing if Sum_imgs( nvox/nv0 ) = numTransformParameters
  unsigned int nVoxPerNVDim = 0;
  for (size_t nImg = 0; nImg < parameterImagesIn.size(); ++nImg)
  {
    nVoxPerNVDim += parameterImagesIn[nImg]->nvox / parameterImagesIn[nImg]->nv;
  }

  if (nVoxPerNVDim != this->transform->GetNumberOfParameters())
  {
    supremo_print_error( "Number of parameters per transformation not as expected." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  // Final sanity check: Total number of voxels equal to number of correspondence model parameters
  unsigned int nTotalVox = 0;
  for (size_t nImg = 0; nImg < parameterImagesIn.size(); ++nImg)
  {
    nTotalVox += parameterImagesIn[nImg]->nvox;
  }

  if (nTotalVox != this->numberOfModelParameters)
  {
    supremo_print_error( "Number of model parameters not as expected." );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  
  // Copy the voxel contents into the parameters
  // Continuously filling parameter array
  unsigned int copyOffsetParameters = 0;
  
  for (unsigned int iSur = 0; iSur < this->numberOfSurrogateSignals; ++iSur)
  {
    for (unsigned int iImg = 0; iImg < parameterImagesIn.size(); ++iImg)
    {
      // Offset from one dimension (v) to the next
      unsigned int dimVSize = parameterImagesIn[iImg]->nvox / parameterImagesIn[iImg]->nv;
      
      // Pointer to the start of the image data 
      PrecisionType* sourcePointer = (PrecisionType*)parameterImagesIn[iImg]->data;
      // forard this to the current v-dimension(curSur ( number of voxels per v-dimension)
      sourcePointer = &(sourcePointer[iSur * dimVSize]);
      
      // Pointer to the model parameters
      PrecisionType* destinationPointer = &(this->modelParameters[copyOffsetParameters]);
      memcpy(destinationPointer, sourcePointer, sizeof( PrecisionType )* dimVSize );
      copyOffsetParameters += dimVSize;
    }
  }
  
  return;
}




//--------------------------------------------
// CorrespondenceModel::GetMaxParameterLength
//--------------------------------------------
CorrespondenceModel::PrecisionType CorrespondenceModel::GetMaxParameterLength( PrecisionType * parametersIn )
{
  PrecisionType maxParameterLength = 0;

  for (unsigned int iSurrSig = 0; iSurrSig < this->numberOfSurrogateSignals; ++iSurrSig)
  {
    PrecisionType* ptrFirstParam = &(parametersIn[iSurrSig * this->numberOfTransformationParameters]);
    PrecisionType curParameterLength = this->transform->GetMaxTransformationParameterLength(ptrFirstParam);
    maxParameterLength = (curParameterLength > maxParameterLength) ? curParameterLength : maxParameterLength;
  }

  return maxParameterLength;
}




//---------------------------------------------------------------------------
// CorrespondenceModel::GetTransformationParameterGradientWRTModelParameters
//---------------------------------------------------------------------------
void CorrespondenceModel::GetTransformationParameterGradientWRTModelParameters( CorrespondenceModel::PrecisionType* transformationGradientIn, 
                                                                                const CorrespondenceModel::SurrogateSignalType & surrogateIn,  
                                                                                CorrespondenceModel::PrecisionType* correspondenceModelGradientOut )
{
  // Need to compute the correspondence model gradient from the transformation parameter gradient
  for (unsigned int iSurr = 0; iSurr < this->numberOfSurrogateSignals; ++iSurr)
  {
    PrecisionType curSurrVal = surrogateIn[iSurr];

    for (unsigned int iTrafo = 0; iTrafo < this->numberOfTransformationParameters; ++iTrafo)
    {
      correspondenceModelGradientOut[iTrafo + iSurr*this->numberOfTransformationParameters] += curSurrVal * transformationGradientIn[iTrafo];
    }
  }
}

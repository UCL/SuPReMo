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
#include "SlidingTransformation.h"




//----------------------------------------------
// SlidingTransformation::SlidingTransformation
//----------------------------------------------
SlidingTransformation::SlidingTransformation( std::shared_ptr<Transformation> transformationAIn, std::shared_ptr<Transformation> transformationBIn,
  nifti_image* defSpaceImage, unsigned int numberOfLevels, unsigned int numberOfLevelsToPerform ) :
  transformA( transformationAIn ),
  transformB( transformationBIn ),
  signedDistMapImage( nullptr ),
  defSpaceImagePyramid( nullptr ),
  lastInitialisedLevel( -1 ),
  signedDistMapImageTransformedA( nullptr ),
  signedDistMapImageTransformedB( nullptr )
{
  // Initialisation of baseclass parameters
  this->interpolation = 1;
  this->deformationVectorFieldImage = nullptr;
  this->warpedPaddingValue = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
  this->dvfImageUpdateRequired = true;

  // Generate the image pyramid that defines in which area the constraint term will be calculated
  this->defSpaceImagePyramid = std::make_shared<ImagePyramid<PrecisionType>>();
  this->defSpaceImagePyramid->GenerateLevels( defSpaceImage, numberOfLevels, numberOfLevelsToPerform );

  // The input transfomrations already know how many parameters they have
  // so need to set this accordingly.
  this->numberOfParameters = this->transformA->GetNumberOfParameters() + this->transformB->GetNumberOfParameters();

  return;
}




//----------------------------------------------
// SlidingTransformation::SlidingTransformation
//----------------------------------------------
SlidingTransformation::SlidingTransformation( const SlidingTransformation& transformToCopy )
{
  // Copy constructor
  // Copy all parameters explicitly

  // Base-class parameters first
  this->dvfImageUpdateRequired = transformToCopy.dvfImageUpdateRequired;
  this->numberOfParameters = transformToCopy.numberOfParameters;
  this->interpolation = transformToCopy.interpolation;
  this->warpedPaddingValue = transformToCopy.warpedPaddingValue;

  // DVF if it exists
  if (transformToCopy.deformationVectorFieldImage != nullptr)
  {
    this->deformationVectorFieldImage = nifti_copy_nim_info( transformToCopy.deformationVectorFieldImage );

    this->deformationVectorFieldImage->data = malloc( deformationVectorFieldImage->nvox * deformationVectorFieldImage->nbyper );
    memcpy( this->deformationVectorFieldImage->data, transformToCopy.deformationVectorFieldImage->data, this->deformationVectorFieldImage->nvox * this->deformationVectorFieldImage->nbyper );
  }
  else
  {
    this->deformationVectorFieldImage = nullptr;
  }

  // Sliding transform specific parfameters

  // Copy the individual transformations first
  this->transformA = transformToCopy.transformA->DeepCopy();
  this->transformB = transformToCopy.transformB->DeepCopy();

  // Copy the def-space image pyramid
  this->defSpaceImagePyramid = transformToCopy.defSpaceImagePyramid;
  this->lastInitialisedLevel = transformToCopy.lastInitialisedLevel;

  // Set the current def space image to the correct level of the pyramid
  if (lastInitialisedLevel >= 0)
  {
    this->curDefSpaceImage = this->defSpaceImagePyramid->GetLevel( this->lastInitialisedLevel );
  }
  else
  {
    this->curDefSpaceImage = nullptr;
  }

  // Copy the signed distance map image
  if (nullptr == transformToCopy.signedDistMapImage)
  {
    this->signedDistMapImage = nullptr;
  }
  else 
  {
    this->signedDistMapImage = nifti_copy_nim_info( transformToCopy.signedDistMapImage );
    this->signedDistMapImage->data = malloc( this->signedDistMapImage->nvox * this->signedDistMapImage->nbyper );
    memcpy( this->signedDistMapImage->data, transformToCopy.signedDistMapImage->data, this->signedDistMapImage->nvox * this->signedDistMapImage->nbyper );
  }
  
  // The transfomred signed distance maps
  if (nullptr != transformToCopy.signedDistMapImageTransformedA)
  {
    this->signedDistMapImageTransformedA = nifti_copy_nim_info( transformToCopy.signedDistMapImageTransformedA );
    this->signedDistMapImageTransformedA->data = malloc( this->signedDistMapImageTransformedA ->nvox * this->signedDistMapImageTransformedA->nbyper );
    memcpy( this->signedDistMapImageTransformedA->data, 
      transformToCopy.signedDistMapImageTransformedA->data, 
      this->signedDistMapImageTransformedA->nvox * this->signedDistMapImageTransformedA->nbyper );
  }
  else
  {
    this->signedDistMapImageTransformedA = nullptr;
  }

  if (nullptr != transformToCopy.signedDistMapImageTransformedB)
  {
    this->signedDistMapImageTransformedB = nifti_copy_nim_info( transformToCopy.signedDistMapImageTransformedB );
    this->signedDistMapImageTransformedB->data = malloc( this->signedDistMapImageTransformedB->nvox * this->signedDistMapImageTransformedB->nbyper );
    memcpy( this->signedDistMapImageTransformedB->data,
      transformToCopy.signedDistMapImageTransformedB->data,
      this->signedDistMapImageTransformedB->nvox * this->signedDistMapImageTransformedB->nbyper );
  }
  else
  {
    this->signedDistMapImageTransformedB = nullptr;
  }

  // The constraint weight value
  this->gapOverlapConstraintWeight = transformToCopy.gapOverlapConstraintWeight;
}




//-----------------------------------------------
// SlidingTransformation::~SlidingTransformation
//-----------------------------------------------
SlidingTransformation::~SlidingTransformation()
{
  // Remove the transfomrs, ...
  this->transformA = nullptr;
  this->transformB = nullptr;
  
  // ... the signed distance maps, ....
  this->ClearSignedDistanceMap();
  this->ClearTransformedSignedDistanceMaps();
  
  // ... and the DVF
  this->ClearDeformationVectorFieldImage();
}




//--------------------------------------------
// SlidingTransformation::GetCopyOfParameters
//--------------------------------------------
SlidingTransformation::PrecisionType* SlidingTransformation::GetCopyOfParameters()
{
  this->numberOfParameters;
  
  // Allocate memory large enough to keep all parameters
  PrecisionType* combinedTransformParameters = (PrecisionType*)malloc( this->GetNumberOfParameters() * sizeof( PrecisionType ) );
  
  // First half of the memory is copied into the first bit
  PrecisionType* parametersTrafoA = this->transformA->GetCopyOfParameters();
  PrecisionType* parametersTrafoB = this->transformB->GetCopyOfParameters();

  // Simple concatenation of parameters
  memcpy( &(combinedTransformParameters[0]), parametersTrafoA, sizeof( PrecisionType ) * this->transformA->GetNumberOfParameters() );
  memcpy( &(combinedTransformParameters[this->transformA->GetNumberOfParameters()]), parametersTrafoB, sizeof( PrecisionType ) * this->transformB->GetNumberOfParameters() );

  // Free the copied memory from individual transformations
  free( parametersTrafoA );
  free( parametersTrafoB );

  return combinedTransformParameters;
}





//----------------------------------------
// SlidingTransformation::InitialiseLevel
//----------------------------------------
void SlidingTransformation::InitialiseLevel( unsigned int levelIn )
{
  // Need to pass on this initialisation to the individual transformatiobns
  this->transformA->InitialiseLevel( levelIn );
  this->transformB->InitialiseLevel( levelIn );
  
  // Handle any sliding transformation specific level-initialisations here
  // Update the current def-space image
  this->curDefSpaceImage = this->defSpaceImagePyramid->GetLevel( levelIn );
  
  // Update the total number of parameters
  this->numberOfParameters = this->transformA->GetNumberOfParameters() + this->transformB->GetNumberOfParameters();

  // Keep a record of which level was initialised
  this->lastInitialisedLevel = levelIn;

  // On each level, an update of the DVF is required.
  this->dvfImageUpdateRequired = true;
}




//--------------------------------------------------
// SlidingTransformation::GetDeformationVectorField
//--------------------------------------------------
nifti_image* SlidingTransformation::GetDeformationVectorField( nifti_image* targetImageIn )
{
  // Only update the calculation if it is required, otherwise return the previously calculated version
  if (false == this->CheckDVFImageUpdateRequired( targetImageIn ))
  {
    return this->deformationVectorFieldImage;
  }
  
  // Calculate a new DVF otherwise
  this->ClearDeformationVectorFieldImage();
  
  // Allocate the transformed signed distance maps, old ones will be cleared 
  this->AllocateTransformedSignedDistanceMaps( targetImageIn );

  // Transform the signed distance maps with both transformations
  this->transformA->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedA );
  this->transformB->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedB );

  // Get the DVFs from the single transfomrs
  nifti_image* dvfTransformA = this->transformA->GetDeformationVectorField( targetImageIn );
  nifti_image* dvfTransformB = this->transformB->GetDeformationVectorField( targetImageIn );

  // Allocate the dvf of the sliding transformation could use A or B...
  this->deformationVectorFieldImage = nifti_copy_nim_info( dvfTransformA );
  this->deformationVectorFieldImage->data = calloc( this->deformationVectorFieldImage->nvox, this->deformationVectorFieldImage->nbyper );

  // Assumption: At this point the DVFs and the warped DVFs have the same number of voxels
  //             except from the different dimensions of course. Potentially need to check this. 

  // 2D implementation 
  // note: decide on dimensionality of DVF, not nz!
  if (dvfTransformA->nu == 2)
  {
    // Pointers to the output DVF
    PrecisionType* dvfPtrX = static_cast<PrecisionType*>(this->deformationVectorFieldImage->data);
    PrecisionType* dvfPtrY = &dvfPtrX[this->deformationVectorFieldImage->nx * this->deformationVectorFieldImage->ny];

    // Pointers to DVF-A (<0)
    PrecisionType* dvfAPtrX = static_cast<PrecisionType*>(dvfTransformA->data);
    PrecisionType* dvfAPtrY = &dvfAPtrX[dvfTransformA->nx * dvfTransformA->ny];

    // Pointers to DVF-B (>=0)
    PrecisionType* dvfBPtrX = static_cast<PrecisionType*>(dvfTransformB->data);
    PrecisionType* dvfBPtrY = &dvfBPtrX[dvfTransformA->nx * dvfTransformA->ny];

   // Pointers to the distance map transformed with A and B
    PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
    PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);

    // Iterate over all voxels
    /// ToDo: make parallel
    for (unsigned int i = 0; i < this->signedDistMapImageTransformedA->nvox; ++i)
    {
      // transformed distance maps (distA/distB) will contain NaN values if the transform
      // maps the voxel outside the extent of the distance map so need to check
      // for NaN values
      if (distAPtr[i] != distAPtr[i])
      {
        if (distBPtr[i] != distBPtr[i])
        {
          // both transformed distance maps are NaN, so the result has to be NaN, too
          dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
        }
        else // dvfB is valid value
        {
          //check if dfvB maps into region A, i.e. if distB < 0
          if (distBPtr[i] < 0)
          {
            // set combined def field to NaN
            dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          }
          else
          {
            // set combined def field to dvfB 
            dvfPtrX[i] = dvfBPtrX[i];
            dvfPtrY[i] = dvfBPtrY[i];
          }
        }
      }
      else 
      {
        // distA is not NaN, but still need to check distB
        if (distBPtr[i] != distBPtr[i])
        {
          //distB is NaN so check if dvfA maps into region B, i.e. if distA >= 0
          if (distAPtr[i] >= 0)
          {
            // set combined def field to NaN
            dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          }
          else 
          {
            // distA < 0 so set combined def field to dvfA
            dvfPtrX[i] = dvfAPtrX[i];
            dvfPtrY[i] = dvfAPtrY[i];
          }
        }
        else 
        {
          // distA and distB are both not NaN
          if (distAPtr[i] + distBPtr[i] < 0)
          {
            // if sum of distMaps < 0 set combined dvf to dvfA 
            dvfPtrX[i] = dvfAPtrX[i];
            dvfPtrY[i] = dvfAPtrY[i];
          }
          else
          {
            // if sum of distMaps >= 0 set combined dvf to dvfB
            dvfPtrX[i] = dvfBPtrX[i];
            dvfPtrY[i] = dvfBPtrY[i];
          }
        }
      }
    }
  }
  else 
  {
    // 3D implementation

    // Pointers to the output DVF
    PrecisionType* dvfPtrX = static_cast<PrecisionType*>(this->deformationVectorFieldImage->data);
    PrecisionType* dvfPtrY = &dvfPtrX[this->deformationVectorFieldImage->nx * this->deformationVectorFieldImage->ny * this->deformationVectorFieldImage->nz];
    PrecisionType* dvfPtrZ = &dvfPtrY[this->deformationVectorFieldImage->nx * this->deformationVectorFieldImage->ny * this->deformationVectorFieldImage->nz];

    // Pointers to the DVF-A
    PrecisionType* dvfAPtrX = static_cast<PrecisionType*>(dvfTransformA->data);
    PrecisionType* dvfAPtrY = &dvfAPtrX[dvfTransformA->nx * dvfTransformA->ny * dvfTransformA->nz];
    PrecisionType* dvfAPtrZ = &dvfAPtrY[dvfTransformA->nx * dvfTransformA->ny * dvfTransformA->nz];

    // Pointers to the DVF-B
    PrecisionType* dvfBPtrX = static_cast<PrecisionType*>(dvfTransformB->data);
    PrecisionType* dvfBPtrY = &dvfBPtrX[dvfTransformB->nx * dvfTransformB->ny * dvfTransformB->nz];
    PrecisionType* dvfBPtrZ = &dvfBPtrY[dvfTransformB->nx * dvfTransformB->ny * dvfTransformB->nz];

    // Pointers to the distance map transformed with A
    PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
    PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);


    // Iterate over all voxels
    // todo: make parallel
    for (unsigned int i = 0; i < this->signedDistMapImageTransformedA->nvox; ++i)
    {
      // transformed distance maps (distA/distB) will contain NaN values if the transform
      // maps the voxel outside the extent of the distance map so need to check
      // for NaN values
      if (distAPtr[i] != distAPtr[i])
      {
        if (distBPtr[i] != distBPtr[i])
        {
          // both transformed distance maps are NaN, so the result has to be NaN, too
          dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          dvfPtrZ[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
        }
        else // dvfB is valid value
        {
          //check if dfvB maps into region A, i.e. if distB < 0
          if (distBPtr[i] < 0)
          {
            // set combined def field to NaN
            dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrZ[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          }
          else
          {
            // set combined def field to dvfB 
            dvfPtrX[i] = dvfBPtrX[i];
            dvfPtrY[i] = dvfBPtrY[i];
            dvfPtrZ[i] = dvfBPtrZ[i];
          }
        }
      }
      else
      {
        // distA is not NaN, but still need to check distB
        if (distBPtr[i] != distBPtr[i])
        {
          //distB is NaN so check if dvfA maps into region B, i.e. if distA >= 0
          if (distAPtr[i] >= 0)
          {
            // set combined def field to NaN
            dvfPtrX[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrY[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
            dvfPtrZ[i] = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
          }
          else
          {
            // distA < 0 so set combined def field to dvfA
            dvfPtrX[i] = dvfAPtrX[i];
            dvfPtrY[i] = dvfAPtrY[i];
            dvfPtrZ[i] = dvfAPtrZ[i];
          }
        }
        else
        {
          // distA and distB are both not NaN
          if (distAPtr[i] + distBPtr[i] < 0)
          {
            // if sum of distMaps < 0 set combined dvf to dvfA 
            dvfPtrX[i] = dvfAPtrX[i];
            dvfPtrY[i] = dvfAPtrY[i];
            dvfPtrZ[i] = dvfAPtrZ[i];
          }
          else
          {
            // if sum of distMaps >= 0 set combined dvf to dvfB
            dvfPtrX[i] = dvfBPtrX[i];
            dvfPtrY[i] = dvfBPtrY[i];
            dvfPtrZ[i] = dvfBPtrZ[i];
          }
        }
      }
    }
  }
  
  // No further update required until parameters or target image etc. were changed
  this->dvfImageUpdateRequired = false;

  return this->deformationVectorFieldImage;
}




//-------------------------------------------------------------------------
// SlidingTransformation::GetConstraintGradientWRTTransformationParameters
//-------------------------------------------------------------------------
SlidingTransformation::PrecisionType* SlidingTransformation::GetConstraintGradientWRTTransformationParameters()
{
  // Allocate the output result
  PrecisionType* outGradWRTTransformParams = (PrecisionType*)malloc( this->numberOfParameters * sizeof( PrecisionType ) );

  // Get the numbers of parameters for both transformations
  size_t numParamsA = this->transformA->GetNumberOfParameters();
  size_t numParamsB = this->transformB->GetNumberOfParameters();
  
  // Compute the constraint gradient of the individual transformations
  PrecisionType* constraintGradientWRTParamsA = this->transformA->GetConstraintGradientWRTTransformationParameters();
  PrecisionType* constraintGradientWRTParamsB = this->transformB->GetConstraintGradientWRTTransformationParameters();


  // Only compute the gap/overlap constraint realted values if the weight is leq 0
  if (this->gapOverlapConstraintWeight > 0)
  {
    // Gap/overlap constraint calculation on a  voxel level first, then convert to gradient wet transformation parameters depending on the individual transformations

    // Allocate a gap/overlap constraint gradient image in line with the DVF size. 
    // Request the current full FOV
    nifti_image* gapOverlapConstraintGradientWRTDVFA = nifti_copy_nim_info( this->transformA->GetDeformationVectorField( this->curDefSpaceImage ) );
    nifti_image* gapOverlapConstraintGradientWRTDVFB = nifti_copy_nim_info( gapOverlapConstraintGradientWRTDVFA );

    // Initialise the DVFS with zero
    gapOverlapConstraintGradientWRTDVFA->data = calloc( gapOverlapConstraintGradientWRTDVFA->nvox, gapOverlapConstraintGradientWRTDVFA->nbyper );
    gapOverlapConstraintGradientWRTDVFB->data = calloc( gapOverlapConstraintGradientWRTDVFB->nvox, gapOverlapConstraintGradientWRTDVFB->nbyper );

    // Calculate the gradients: 
    // Note, the function GetImageGradientWRTDVF already reorientates the image gradient into real-world space
    // so no reorientation required below.
    // dD_A / dT_A 
    this->transformA->GetImageGradientWRTDVF( this->signedDistMapImage, gapOverlapConstraintGradientWRTDVFA );
    // dD_B / dT_B 
    this->transformB->GetImageGradientWRTDVF( this->signedDistMapImage, gapOverlapConstraintGradientWRTDVFB );

    // Check if the previously calculated signed distance maps can be used
    // To make the constraint term independent of the target image size, here the 
    // signed distance map is used. 
    if (!this->TransformedSignedDistMapsCanBeUsedForConstraintCalculations())
    {
      // Allocate the the transformed signed distance map
      this->AllocateTransformedSignedDistanceMaps( this->curDefSpaceImage );

      // Use the individual transformations to warp them 
      this->transformA->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedA );
      this->transformB->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedB );
    }

    // Calculate 
    // dC / dT_A = - D_B * (dD_A / dT_A)   if: D_A * D_B < 0 
    // dC / dT_B = - D_A * (dD_B / dT_B)   if: D_A * D_B < 0 

    if (gapOverlapConstraintGradientWRTDVFA->nu == 2)
    {
      // 2D implementation
      size_t numVoxels = this->signedDistMapImageTransformedA->nx * this->signedDistMapImageTransformedA->ny;

      // Pointers to the signed distance maps D_A and D_B  
      PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
      PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);

      // Pointers to the gap-overlap constraint gradients A and B; X and Y
      PrecisionType* goctGradAPtrX = static_cast<PrecisionType*>(gapOverlapConstraintGradientWRTDVFA->data);
      PrecisionType* goctGradAPtrY = &goctGradAPtrX[numVoxels];

      PrecisionType* goctGradBPtrX = static_cast<PrecisionType*>(gapOverlapConstraintGradientWRTDVFB->data);
      PrecisionType* goctGradBPtrY = &goctGradBPtrX[numVoxels];

      // iterate over voxels
      for (size_t i = 0; i < numVoxels; ++i)
      {
        if (distAPtr[i] * distBPtr[i] < 0)
        {
          // Multiply the distance map gradient of transform A and B with -D_B or -D_A respectively
          goctGradAPtrX[i] *= -distBPtr[i];
          goctGradAPtrY[i] *= -distBPtr[i];

          goctGradBPtrX[i] *= -distAPtr[i];
          goctGradBPtrY[i] *= -distAPtr[i];
        }
        else
        {
          // Set the constraint gradient to zero
          goctGradAPtrX[i] = static_cast<PrecisionType>(0.);
          goctGradAPtrY[i] = static_cast<PrecisionType>(0.);

          goctGradBPtrX[i] = static_cast<PrecisionType>(0.);
          goctGradBPtrY[i] = static_cast<PrecisionType>(0.);
        }
      }
    }
    else
    {
      // 3D implementation

      size_t numVoxels = this->signedDistMapImageTransformedA->nx * this->signedDistMapImageTransformedA->ny * this->signedDistMapImageTransformedA->nz;

      // Pointers to the signed distance maps D_A and D_B  
      PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
      PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);

      // Pointers to the gap-overlap constraint gradients A and B; X and Y
      PrecisionType* goctGradAPtrX = static_cast<PrecisionType*>(gapOverlapConstraintGradientWRTDVFA->data);
      PrecisionType* goctGradAPtrY = &goctGradAPtrX[numVoxels];
      PrecisionType* goctGradAPtrZ = &goctGradAPtrY[numVoxels];

      PrecisionType* goctGradBPtrX = static_cast<PrecisionType*>(gapOverlapConstraintGradientWRTDVFB->data);
      PrecisionType* goctGradBPtrY = &goctGradBPtrX[numVoxels];
      PrecisionType* goctGradBPtrZ = &goctGradBPtrY[numVoxels];

      // iterate over voxels
      for (size_t i = 0; i < numVoxels; ++i)
      {
        if (distAPtr[i] * distBPtr[i] < 0)
        {
          // Multiply the distance map gradient of transform A and B with -D_B or -D_A respectively
          goctGradAPtrX[i] *= -distBPtr[i];
          goctGradAPtrY[i] *= -distBPtr[i];
          goctGradAPtrZ[i] *= -distBPtr[i];

          goctGradBPtrX[i] *= -distAPtr[i];
          goctGradBPtrY[i] *= -distAPtr[i];
          goctGradBPtrZ[i] *= -distAPtr[i];
        }
        else
        {
          // Set the constraint gradient to zero
          goctGradAPtrX[i] = static_cast<PrecisionType>(0.);
          goctGradAPtrY[i] = static_cast<PrecisionType>(0.);
          goctGradAPtrZ[i] = static_cast<PrecisionType>(0.);

          goctGradBPtrX[i] = static_cast<PrecisionType>(0.);
          goctGradBPtrY[i] = static_cast<PrecisionType>(0.);
          goctGradBPtrZ[i] = static_cast<PrecisionType>(0.);
        }
      }
    }

    // Convert the gradient to transformation parameters
    PrecisionType* gapOverlapConstraintGradientWRTParamsTA = this->transformA->GetDVFGradientWRTTransformationParameters( gapOverlapConstraintGradientWRTDVFA );
    PrecisionType* gapOverlapConstraintGradientWRTParamsTB = this->transformB->GetDVFGradientWRTTransformationParameters( gapOverlapConstraintGradientWRTDVFB );

    // First half of parameters get results of A
    for (size_t i = 0; i < numParamsA; ++i)
    {
      outGradWRTTransformParams[i] = constraintGradientWRTParamsA[i] + this->gapOverlapConstraintWeight * gapOverlapConstraintGradientWRTParamsTA[i];
    }
    // Second half gets results of B
    for (size_t i = 0; i < numParamsB; ++i)
    {
      outGradWRTTransformParams[i + numParamsA] = constraintGradientWRTParamsB[i] + this->gapOverlapConstraintWeight * gapOverlapConstraintGradientWRTParamsTB[i];
    }
    
    // Free the images that were only used within this scope
    nifti_image_free( gapOverlapConstraintGradientWRTDVFA );
    nifti_image_free( gapOverlapConstraintGradientWRTDVFB );

    free( gapOverlapConstraintGradientWRTParamsTA );
    free( gapOverlapConstraintGradientWRTParamsTB );
  }
  else // this->gapOverlapConstraintWeight > 0, so we can ignore the gap/overlap contriubtion to the gradient
  {
    // First half of parameters get results of A
    for (size_t i = 0; i < numParamsA; ++i)
    {
      outGradWRTTransformParams[i] = constraintGradientWRTParamsA[i];
    }
    // Second half gets results of B
    for (size_t i = 0; i < numParamsB; ++i)
    {
      outGradWRTTransformParams[i + numParamsA] = constraintGradientWRTParamsB[i];
    }
  }

  // Free the individual results
  free( constraintGradientWRTParamsA );
  free( constraintGradientWRTParamsB );

  return outGradWRTTransformParams;
}




//-------------------------------------------
// SlidingTransformation::GetConstraintValue
//-------------------------------------------
double SlidingTransformation::GetConstraintValue()
{
#ifdef _DEBUG
  std::cout << "Called SlidingTransformation::GetConstraintValue()" << std::endl;
#endif
  // Get the constraint values for each transformation individually
  double constraintA = this->transformA->GetConstraintValue();
  double constraintB = this->transformA->GetConstraintValue();

  // Get the sliding-specific constraint value
  double gapOverlapConstraint = this->CalculateGapOverlapConstraintTerm();
 
  // then sum everything up and return this.
  double constraintSum = constraintA + constraintB + gapOverlapConstraint;

  return constraintSum;
}




//------------------------------------------------------------------
// SlidingTransformation::GetDVFGradientWRTTransformationParameters
//------------------------------------------------------------------
SlidingTransformation::PrecisionType* SlidingTransformation::GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn )
{
  // Need to split the input dense DVF in into the two regions, 
  // then get the gradient from the individual transformations.
  

  // need to check if we can still use the deformed signed distance maps, otherwise need to update 
  if (!this->CheckImageGeometryEquality( this->signedDistMapImageTransformedA, denseDVFIn ))
  {
    // Update the signed distance maps
    nifti_image* scalarValuedImageFromDVF = nifti_copy_nim_info( denseDVFIn );
    
    // Make sure the image is scalar valued
    int vectorDim = scalarValuedImageFromDVF->nu;
    scalarValuedImageFromDVF->nu = 1;
    scalarValuedImageFromDVF->dim[5] = 1;
    scalarValuedImageFromDVF->dim[0] = vectorDim;
    scalarValuedImageFromDVF->ndim = vectorDim;
    scalarValuedImageFromDVF->nvox = scalarValuedImageFromDVF->nx * scalarValuedImageFromDVF->ny * scalarValuedImageFromDVF->nz;

    // Generate the transformed signed distance maps
    this->AllocateTransformedSignedDistanceMaps( scalarValuedImageFromDVF );
    this->transformA->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedA );
    this->transformB->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedB );
    
    nifti_image_free( scalarValuedImageFromDVF );
  }
   
  // Split the DVF
  nifti_image* denseDVFA = nifti_copy_nim_info( denseDVFIn );
  nifti_image* denseDVFB = nifti_copy_nim_info( denseDVFIn );

  // initialise with zero
  denseDVFA->data = calloc( denseDVFA->nvox, denseDVFA->nbyper );
  denseDVFB->data = calloc( denseDVFB->nvox, denseDVFB->nbyper );

  // note: decide on dimensionality of DVF->nu, not nz!
  if (denseDVFIn->nu == 2)
  {
    // 2D implementation
    size_t numVox = denseDVFIn->nx * denseDVFIn->ny;

    // Pointers to the distance map transformed with A and B
    PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
    PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);

    // Pointers to dense DVFs A and B
    PrecisionType* dvfAPtrX = static_cast<PrecisionType*>(denseDVFA->data);
    PrecisionType* dvfAPtrY = &dvfAPtrX[numVox];

    PrecisionType* dvfBPtrX = static_cast<PrecisionType*>(denseDVFB->data);
    PrecisionType* dvfBPtrY = &dvfBPtrX[numVox];

    // Pointers to dense (combined/complete) DVFs
    PrecisionType* dvfInPtrX = static_cast<PrecisionType*>(denseDVFIn->data);
    PrecisionType* dvfInPtrY = &dvfInPtrX[numVox];

    // iterate over all voxels
    for (unsigned int i = 0; i < numVox; ++i)
    {
      // need to check for NaNs in distance maps D_A and D_B
      // if D_A is NaN and D_B >= 0 (indicating region-B transform maps into region-B)
      // then copy voxel-based gradient value in to denseDVFB
      if (distAPtr[i] != distAPtr[i] && distBPtr[i] >= 0)
      {
        dvfBPtrX[i] = dvfInPtrX[i];
        dvfBPtrY[i] = dvfInPtrY[i];
      }

      // if D_B is NaN and D_A < 0 (indicating region-A transform maps into region-A)
      // then copy voxel-based gradient value in to A
      if (distBPtr[i] != distBPtr[i] && distAPtr[i] < 0)
      {
        dvfAPtrX[i] = dvfInPtrX[i];
        dvfAPtrY[i] = dvfInPtrY[i];
      }

      // if both SDMs (D_A and D_B) are not NaN then assign voxel-based gradient value to correct region
      // based on D_A + D_B, otherwise leave at zero to have no contribution. 
      if (distAPtr[i] == distAPtr[i] && distBPtr[i] == distBPtr[i])
      {
        if (distAPtr[i] + distBPtr[i] < 0)
        {
          dvfAPtrX[i] = dvfInPtrX[i];
          dvfAPtrY[i] = dvfInPtrY[i];
        }
        else 
        {
          dvfBPtrX[i] = dvfInPtrX[i];
          dvfBPtrY[i] = dvfInPtrY[i];
        }
      }
    }
  }
  else 
  {
    // 3D implementation
    size_t numVox = denseDVFIn->nx * denseDVFIn->ny * denseDVFIn->nz;

    // Pointers to the distance map transformed with A and B
    PrecisionType* distAPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedA->data);
    PrecisionType* distBPtr = static_cast<PrecisionType*>(this->signedDistMapImageTransformedB->data);

    // Pointers to dense DVFs A and B
    PrecisionType* dvfAPtrX = static_cast<PrecisionType*>(denseDVFA->data);
    PrecisionType* dvfAPtrY = &dvfAPtrX[numVox];
    PrecisionType* dvfAPtrZ = &dvfAPtrY[numVox];

    PrecisionType* dvfBPtrX = static_cast<PrecisionType*>(denseDVFB->data);
    PrecisionType* dvfBPtrY = &dvfBPtrX[numVox];
    PrecisionType* dvfBPtrZ = &dvfBPtrY[numVox];

    // Pointers to dense (combined/complete) DVFs
    PrecisionType* dvfInPtrX = static_cast<PrecisionType*>(denseDVFIn->data);
    PrecisionType* dvfInPtrY = &dvfInPtrX[numVox];
    PrecisionType* dvfInPtrZ = &dvfInPtrY[numVox];

    // iterate over all voxels
    for (unsigned int i = 0; i < numVox; ++i)
    {
      // need to check for NaNs in distance maps D_A and D_B
      // if D_A is NaN and D_B >= 0 (indicating region-B transform maps into region-B)
      // then copy voxel-based gradient value in to denseDVFB
      if (distAPtr[i] != distAPtr[i] && distBPtr[i] >= 0)
      {
        dvfBPtrX[i] = dvfInPtrX[i];
        dvfBPtrY[i] = dvfInPtrY[i];
        dvfBPtrZ[i] = dvfInPtrZ[i];
      }

      // if D_B is NaN and D_A < 0 (indicating region-A transform maps into region-A)
      // then copy voxel-based gradient value in to A
      if (distBPtr[i] != distBPtr[i] && distAPtr[i] < 0)
      {
        dvfAPtrX[i] = dvfInPtrX[i];
        dvfAPtrY[i] = dvfInPtrY[i];
        dvfAPtrZ[i] = dvfInPtrZ[i];
      }

      // if both SDMs (D_A and D_B) are not NaN then assign voxel-based gradient value to correct region
      // based on D_A + D_B, otherwise leave at zero to have no contribution. 
      if (distAPtr[i] == distAPtr[i] && distBPtr[i] == distBPtr[i])
      {
        if (distAPtr[i] + distBPtr[i] < 0)
        {
          dvfAPtrX[i] = dvfInPtrX[i];
          dvfAPtrY[i] = dvfInPtrY[i];
          dvfAPtrZ[i] = dvfInPtrZ[i];
        }
        else
        {
          dvfBPtrX[i] = dvfInPtrX[i];
          dvfBPtrY[i] = dvfInPtrY[i];
          dvfBPtrZ[i] = dvfInPtrZ[i];
        }
      }
    }
  }

  // Get the gradient from the individual transformations
  PrecisionType* dvfGradWRTParamsA = this->transformA->GetDVFGradientWRTTransformationParameters( denseDVFA );
  PrecisionType* dvfGradWRTParamsB = this->transformB->GetDVFGradientWRTTransformationParameters( denseDVFB );
  
  // Concatenate the result
  PrecisionType* concatenatedResult = static_cast<PrecisionType*>(calloc( this->numberOfParameters, sizeof( PrecisionType ) ));
  
  memcpy( &(concatenatedResult[0]), dvfGradWRTParamsA, sizeof( PrecisionType )* this->transformA->GetNumberOfParameters() );
  memcpy( &(concatenatedResult[this->transformA->GetNumberOfParameters()]), dvfGradWRTParamsB, sizeof( PrecisionType )* this->transformB->GetNumberOfParameters() );
  
  // Free up the result from the individual transforms
  free( dvfGradWRTParamsA );
  free( dvfGradWRTParamsB );
  
  nifti_image_free( denseDVFA );
  nifti_image_free( denseDVFB );

  return concatenatedResult;
}




//------------------------------------------------------
// SlidingTransformation::SetGapOverlapConstraintWeight
//------------------------------------------------------
void SlidingTransformation::SetGapOverlapConstraintWeight( double gapOverlapConstraintIn)
{
  if (gapOverlapConstraintIn < 0. )
  {
    char msg[200];
    sprintf_s( msg, "Gap-overlap constraint weight has to be positive." );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }

  this->gapOverlapConstraintWeight = gapOverlapConstraintIn;
  return;
}




//--------------------------------------
// SlidingTransformation::SetParameters
//--------------------------------------
void SlidingTransformation::SetParameters( SlidingTransformation::PrecisionType* parametersIn, bool parametersAreDisplacements )
{
  // Pass setting of the parameters to the individual transformations
  // transformA --> no offset required
  SlidingTransformation::PrecisionType* pointerToParamsA = parametersIn;

  // transformB --> offset by number of parameters of transformA
  SlidingTransformation::PrecisionType* pointerToParamsB = &(parametersIn[this->transformA->GetNumberOfParameters()]);
  
  this->transformA->SetParameters( pointerToParamsA, parametersAreDisplacements );
  this->transformB->SetParameters( pointerToParamsB, parametersAreDisplacements );
  
  // Force re-calculation of DVF with new parameters
  this->dvfImageUpdateRequired = true;
}





//-----------------------------------------------
// SlidingTransformation::SetBendingEnergyWeight
//-----------------------------------------------
SlidingTransformation::PrecisionType SlidingTransformation::GetSumOfPenaltyWeights()
{
  PrecisionType sumOfPenaltyWeightsA = this->transformA->GetSumOfPenaltyWeights();
  PrecisionType sumOfPenaltyWeightsB = this->transformB->GetSumOfPenaltyWeights();
  
  if (sumOfPenaltyWeightsA != sumOfPenaltyWeightsB)
  {
    supremo_print_warning( "Constraint weights of transformations for region A and B are different. Using average for further computations" );
  }

  PrecisionType sumOfPenaltyWeights = this->gapOverlapConstraintWeight + (sumOfPenaltyWeightsA + sumOfPenaltyWeightsB) / 2.f;
  return sumOfPenaltyWeights;
}




//------------------------------------------------------------
// SlidingTransformation::GetMaxTransformationParameterLength
//------------------------------------------------------------
SlidingTransformation::PrecisionType SlidingTransformation::GetMaxTransformationParameterLength( PrecisionType * parametersIn )
{
  // Pass the parameters to the individual transformations
  // transformA --> no offset required
  SlidingTransformation::PrecisionType* pointerToParamsA = parametersIn;

  // transformB --> offset by number of parameters of transformA
  SlidingTransformation::PrecisionType* pointerToParamsB = &(parametersIn[this->transformA->GetNumberOfParameters()]);

  PrecisionType maxTransformParameterLengthA = this->transformA->GetMaxTransformationParameterLength( pointerToParamsA );
  PrecisionType maxTransformParameterLengthB = this->transformB->GetMaxTransformationParameterLength( pointerToParamsB );

  // Return the maximum of both values
  return maxTransformParameterLengthA > maxTransformParameterLengthB ? maxTransformParameterLengthA : maxTransformParameterLengthB;
}





//-------------------------------------------------
// SlidingTransformation::GetTransformationAsImage
//-------------------------------------------------
std::vector<nifti_image*> SlidingTransformation::GetTransformationAsImage()
{
  // The sliding transformation collects all transformation images from transformation A, and then from B
  // Images are returned in this order A1, A2, A3, ..., AN, B1, B2, ..., BM
  
  std::vector<nifti_image*> retTransformImages;

  std::vector<nifti_image*> trafoAImgs = this->transformA->GetTransformationAsImage();
  std::vector<nifti_image*> trafoBImgs = this->transformB->GetTransformationAsImage();

  for (size_t i = 0; i < trafoAImgs.size(); ++i)
  {
    retTransformImages.push_back( trafoAImgs[i] );
  }

  for (size_t i = 0; i < trafoBImgs.size(); ++i)
  {
    retTransformImages.push_back( trafoBImgs[i] );
  }

  return retTransformImages;
}




//---------------------------------
// SlidingTransformation::DeepCopy
//---------------------------------
std::shared_ptr<Transformation> SlidingTransformation::DeepCopy()
{
  std::shared_ptr<SlidingTransformation> copiedTransform = std::make_shared<SlidingTransformation>( *this );
  return copiedTransform;
}




//--------------------------------------------------------
// SlidingTransformation::DisplayTransformationParameters
//--------------------------------------------------------
void SlidingTransformation::DisplayTransformationParameters()
{
  printf( "Sliding-transformation parameters\n" );
  return;
}




//---------------------------------------------
// SlidingTransformation::SetSignedDistanceMap
//---------------------------------------------
void SlidingTransformation::SetSignedDistanceMapImage( nifti_image* signedDistanceMapImageIn )
{
  this->ClearSignedDistanceMap();

  // Set the signed distance map pyramid to the member variable
  this->signedDistMapImage = nifti_copy_nim_info( signedDistanceMapImageIn );
  this->signedDistMapImage->data = malloc( this->signedDistMapImage->nvox * this->signedDistMapImage->nbyper );
  memcpy( this->signedDistMapImage->data, signedDistanceMapImageIn->data, this->signedDistMapImage->nvox * this->signedDistMapImage->nbyper );

  // The signed distance map will affect how the dvf is calculated. So need to force the update.
  this->dvfImageUpdateRequired = true;
}




//-----------------------------------------------
// SlidingTransformation::ClearSignedDistanceMap
//-----------------------------------------------
void SlidingTransformation::ClearSignedDistanceMap()
{
  // Free the copied nifti image
  if (nullptr != this->signedDistMapImage)
  {
    nifti_image_free( this->signedDistMapImage );
    this->signedDistMapImage = nullptr;
  }
}




//-----------------------------------------------------------
// SlidingTransformation::ClearTransformedSignedDistanceMaps
//-----------------------------------------------------------
void SlidingTransformation::ClearTransformedSignedDistanceMaps()
{
  if (nullptr != this->signedDistMapImageTransformedA)
  {
    nifti_image_free( this->signedDistMapImageTransformedA );
    this->signedDistMapImageTransformedA = nullptr;
  }
  if (nullptr != this->signedDistMapImageTransformedB)
  {
    nifti_image_free( this->signedDistMapImageTransformedB );
    this->signedDistMapImageTransformedB = nullptr;
  }
}




//--------------------------------------------------------------
// SlidingTransformation::AllocateTransformedSignedDistanceMaps
//--------------------------------------------------------------
void SlidingTransformation::AllocateTransformedSignedDistanceMaps( nifti_image* targetImageIn )
{
  // Clean up the member variables before allocating new ones
  this->ClearTransformedSignedDistanceMaps();

  // Transform the signed distance map with transformation A and B
  // Allocate images for the warped SDTs
  this->signedDistMapImageTransformedA = nifti_copy_nim_info( targetImageIn );
  this->signedDistMapImageTransformedB = nifti_copy_nim_info( targetImageIn );

  // make sure the trasnformed dsitance map images are of floating precision type
  this->signedDistMapImageTransformedA->nbyper = sizeof( PrecisionType );
  this->signedDistMapImageTransformedB->nbyper = sizeof( PrecisionType );
  if (sizeof( PrecisionType ) == sizeof( float ))
  {
    this->signedDistMapImageTransformedA->datatype = NIFTI_TYPE_FLOAT32;
    this->signedDistMapImageTransformedB->datatype = NIFTI_TYPE_FLOAT32;
  }
  else
  {
    this->signedDistMapImageTransformedA->datatype = NIFTI_TYPE_FLOAT64;
    this->signedDistMapImageTransformedB->datatype = NIFTI_TYPE_FLOAT64;
  }
  this->signedDistMapImageTransformedA->data = calloc( this->signedDistMapImageTransformedA->nvox, sizeof( PrecisionType ) );
  this->signedDistMapImageTransformedB->data = calloc( this->signedDistMapImageTransformedB->nvox, sizeof( PrecisionType ) );
}




//----------------------------------------------------------
// SlidingTransformation::CalculateGapOverlapConstraintTerm
//----------------------------------------------------------
double SlidingTransformation::CalculateGapOverlapConstraintTerm()
{
  // NOTE - this method assumes the current warped distance maps have already
  // been calculated by calling the GetDeformationField() method prior to calling this
  // method. The GetDeformtionField method will usually be called when warping the image
  // to calculate the image similarities, so this prevents re-calculating the distance maps
  // unnecessarily, but if the image similarities all have a weight of 0 and therefore
  // the warped image is not calculated, the GetDeformationField() method must still be
  // called.

  // NOTE2 - the gap-overlap penalty term is calculated at all voxels within the reference
  // image, even if they are outside the mask or have a NaN value in the reference or
  // warped image - this is to ensure the transformations for the 2 regions are free of
  // gaps and overlaps, even in areas where the images are not being used to drive the
  // registration
  
  if (this->gapOverlapConstraintWeight <= 0)
  {
    return 0.;
  }

  // loop over all voxels and sum up gap-overlap penalty term values from each voxel.
  // the gap-overlap penalty term is defined as -WDM1*WDM2 if WDM1*WDM2<0 (i.e. the
  // WDMs point to different regions, indicating a gap or overlap), and 0 otherwise
  double gapOverlapTotal = 0.;
  double gapOverlapValue = 0.;


  // Check if the previously calculated signed distance maps can be used
  // To make the constraint term independent of the image target size, here the 
  // signed distance map is used. If using full images and a correctly sized 
  // SDM, no recalculation should be required
  if (!this->TransformedSignedDistMapsCanBeUsedForConstraintCalculations())
  {
    // Allocate the the transformed signed distance map
    this->AllocateTransformedSignedDistanceMaps( this->curDefSpaceImage );

    // Use the individual transformations to warp them 
    this->transformA->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedA );
    this->transformB->TransformImage( this->signedDistMapImage, this->signedDistMapImageTransformedB );
  }

  PrecisionType* distAPtr = static_cast<PrecisionType *>(this->signedDistMapImageTransformedA->data);
  PrecisionType *distBPtr = static_cast<PrecisionType *>(this->signedDistMapImageTransformedB->data);

  size_t numVox = this->signedDistMapImageTransformedA->nx * this->signedDistMapImageTransformedA->ny * this->signedDistMapImageTransformedA->nz;

  for (size_t n = 0; n < numVox; ++n)
  {
    gapOverlapValue = distAPtr[n] * distBPtr[n];
    // if NaN value in either WDM then gapOverlapValue = NaN, so will fail
    // test for less than 0
    if (gapOverlapValue < 0)
    {
      gapOverlapTotal -= gapOverlapValue;
    }
  }
  // normalise by the number of voxels and return weighted value
  // gapOverlapTotal /= double(numVox);
  return double( this->gapOverlapConstraintWeight) * gapOverlapTotal;
}








//------------------------------------------------------------------------------------
// SlidingTransformation::TransformedSignedDistMapsCanBeUsedForConstraintCalculations
//------------------------------------------------------------------------------------
bool SlidingTransformation::TransformedSignedDistMapsCanBeUsedForConstraintCalculations()
{
  // Do a simple geometry check between the exisiting transformed singed distance map and the image currently defining the space
  if (nullptr != this->signedDistMapImageTransformedA
    && this->CheckImageGeometryEquality( this->curDefSpaceImage, this->signedDistMapImageTransformedA ))
  {
    return true;
  }
  else
  {
    return false;
  }
}

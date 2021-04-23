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
#include "MathAdditions.h"
#include "Transformation.h"
#include "_reg_resampling.h"
#include "_reg_maths.h"
#include <iostream>




//--------------------------------
// Transformation::TransformImage
//--------------------------------
void Transformation::TransformImage( nifti_image* sourceImgIn, nifti_image* warpedImgIn )
{
  // Calculate the DVF that is the size of the warped image - not more necessary
  this->GetDeformationVectorField( warpedImgIn );
  
  // Perform the resampling according to the freshly calculated DVF
  // Note: The function reg_resampleImage appears suitable, but is not, since
  //       Image dimensionality is determiend by the number of z-slices". For
  //       axial multi-slice acquisitions this test will be insufficient. In 
  //       reg-resp dimensionality was determined by vector dimensionality of
  //       each point: (nu == 3), or (nu == 2). This is replicated below. 
  //       Potentially, at a later stage the functions need to be migrated to
  //       this transformation class. 
  // Note2: If the DTI processing is required, then this has to be implemented. 
  
  
  //reg_resampleImage( sourceImgIn,
  //  warpedImgIn,
  //  dvfImg,
  //  nullptr,
  //  this->interpolation,
  //  this->warpedPaddingValue );

  // cubic interpolation
  if (3 == this->interpolation)
  {
    // 3D
    if (3 == this->deformationVectorFieldImage->nu)
    {
      this->CubicSplineTransformImage3D( sourceImgIn, warpedImgIn );
    }
    else
    {
      // 2D
      this->CubicSplineTransformImage2D( sourceImgIn, warpedImgIn );
    }
  }
  // nearest neighbour interpolation
  else if (0 == this->interpolation)
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->NearestNeighbourTransformImage3D( sourceImgIn, warpedImgIn );
    }
    else
    {
      // 2D
      this->NearestNeighbourTransformImage2D( sourceImgIn, warpedImgIn );
    }
  }
  // tri-/bi-linear interpolation by default
  else
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->TrilinearTransformImage( sourceImgIn, warpedImgIn );
    }
    else
    {
      // 2D
      this->BilinearTransformImage( sourceImgIn, warpedImgIn );
    }
  }
}




//---------------------------------------
// Transformation::TransformImageAdjoint
//---------------------------------------
void Transformation::TransformImageAdjoint( nifti_image * sourceImage, nifti_image * sourceWeightsImage,
  nifti_image * warpedImage, nifti_image * warpedWeightsImage )
{
  // Need to calculate the DVF for each point of the source image
  // note the member variable will be used, so ignore return value
  this->GetDeformationVectorField( sourceImage );

  // need to determine the interpolation method and dimensionality to call correct function
  
  // cubic interpolation
  if (3 == this->interpolation)
  {
    // 3D
    if (3 == this->deformationVectorFieldImage->nu) 
    { 
      this->CubicSplineTransformImageAdjoint3D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else
    {
      // 2D
      this->CubicSplineTransformImageAdjoint2D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
  // nearest neighbour interpolation
  else if (0 == this->interpolation) 
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->NearestNeighbourTransformImageAdjoint3D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else 
    {
      // 2D
      this->NearestNeighbourTransformImageAdjoint2D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
  // tri-/bi-linear interpolation by default
  else 
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->TrilinearTransformImageAdjoint( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else
    {
      // 2D
      this->BilinearTransformImageAdjoint( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
}



//----------------------------------------
// Transformation::GetImageGradientWRTDVF
//----------------------------------------
void Transformation::GetImageGradientWRTDVF( nifti_image* sourceImage, nifti_image* outWarpedGradientImage )
{
  // Calculate the DVF for the size of the output warped gradient image
  this->GetDeformationVectorField( outWarpedGradientImage );
  
  // Calculate the actual gradient
  // Apparently we could have used the below from the nifty-reg library, but 
  // the internal dimensionality switching does not work for multi-slice image
  // acquisitions (including 4dCT and MRI). Hence we need to reimplement here 
  // what as previously available in reg-resp
  //reg_getImageGradient( sourceImage, 
  //                      outWarpedGradientImage, 
  //                      this->deformationVectorFieldImage, 
  //                      nullptr, 
  //                      this->interpolation, 
  //                      this->warpedPaddingValue, 0 );
    // cubic interpolation
  if (3 == this->interpolation)
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->CubicSplineImageGradient3D( sourceImage, outWarpedGradientImage );
    }
    else
    {
      // 2D
      this->CubicSplineImageGradient2D( sourceImage, outWarpedGradientImage );
    }
  }
  // tri-/bi-linear interpolation by default
  else
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->TrilinearImageGradient( sourceImage, outWarpedGradientImage );
    }
    else
    {
      // 2D
      this->BilinearImageGradient( sourceImage, outWarpedGradientImage );
    }
  }

  // Note: This reorientation from index to real-world coordinates was, for performance
  //       reasons, in niftyreg-resp done after calculating the gradient with respect to 
  //       the transformation parameters (less dense = less number of computations). For
  //       clarity and readability, we moved this reorientation here. 
  this->ReorientateVectorImage( outWarpedGradientImage, sourceImage->sto_ijk );
}




//----------------------------------------
// Transformation::TrilinearImageGradient
//----------------------------------------
void Transformation::TrilinearImageGradient(nifti_image* sourceImage, nifti_image* resultGradientImage )
{
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)resultGradientImage->nx*resultGradientImage->ny*resultGradientImage->nz;
  long sourceVoxelNumber = (long)sourceImage->nx * sourceImage->ny * sourceImage->nz;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)resultGradientImage->nx*resultGradientImage->ny*resultGradientImage->nz;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultGradientImagePtr = static_cast<PrecisionType *>(resultGradientImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[targetVoxelNumber];

  PrecisionType paddingValue = this->warpedPaddingValue;

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  // Iteration over the different volume along the 4th axis
  for (int t = 0; t < resultGradientImage->nt; t++)
  {
    PrecisionType *resultGradientPtrX = &resultGradientImagePtr[t * 3 * targetVoxelNumber];
    PrecisionType *resultGradientPtrY = &resultGradientPtrX[targetVoxelNumber];
    PrecisionType *resultGradientPtrZ = &resultGradientPtrY[targetVoxelNumber];

    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    int previous[3], a, b, c, X, Y, Z;
    PrecisionType position[3], xBasis[2], yBasis[2], zBasis[2];
    PrecisionType deriv[2];
    deriv[0] = -1;
    deriv[1] = 1;
    PrecisionType relative, world[3], grad[3], coeff;
    PrecisionType xxTempNewValue, yyTempNewValue, zzTempNewValue, xTempNewValue, yTempNewValue;
    PrecisionType *zPointer, *xyzPointer;
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, world, position, previous, xBasis, yBasis, zBasis, relative, grad, coeff, \
              a, b, c, X, Y, Z, zPointer, xyzPointer, xTempNewValue, yTempNewValue, xxTempNewValue, yyTempNewValue, zzTempNewValue) \
      shared(sourceIntensity, targetVoxelNumber, sourceVoxelNumber, deriv, paddingValue, \
             deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
             sourceIJKMatrix, sourceImage, resultGradientPtrX, resultGradientPtrY, resultGradientPtrZ)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {

      grad[0] = 0.0;
      grad[1] = 0.0;
      grad[2] = 0.0;

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      world[2] = (PrecisionType)deformationFieldPtrZ[index];

      /* real -> voxel; source space */
      reg_mat44_mul( sourceIJKMatrix, world, position );

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      previous[2] = nmm_floorInt( position[2] );
      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;
      // basis values along the z axis
      relative = position[2] - (PrecisionType)previous[2];
      zBasis[0] = (PrecisionType)(1.0 - relative);
      zBasis[1] = relative;

      // The padding value is used for interpolation if it is different from NaN
      if (paddingValue == paddingValue)
      {
        for (c = 0; c < 2; c++)
        {
          Z = previous[2] + c;
          if (Z > -1 && Z < sourceImage->nz)
          {
            zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
            xxTempNewValue = 0.0;
            yyTempNewValue = 0.0;
            zzTempNewValue = 0.0;
            for (b = 0; b < 2; b++)
            {
              Y = previous[1] + b;
              if (Y > -1 && Y < sourceImage->ny)
              {
                xyzPointer = &zPointer[Y*sourceImage->nx + previous[0]];
                xTempNewValue = 0.0;
                yTempNewValue = 0.0;
                for (a = 0; a < 2; a++)
                {
                  X = previous[0] + a;
                  if (X > -1 && X < sourceImage->nx)
                  {
                    coeff = *xyzPointer;
                    xTempNewValue += coeff * deriv[a];
                    yTempNewValue += coeff * xBasis[a];
                  } // end X in range
                  else
                  {
                    xTempNewValue += paddingValue * deriv[a];
                    yTempNewValue += paddingValue * xBasis[a];
                  }
                  xyzPointer++;
                } // end a
                xxTempNewValue += xTempNewValue * yBasis[b];
                yyTempNewValue += yTempNewValue * deriv[b];
                zzTempNewValue += yTempNewValue * yBasis[b];
              } // end Y in range
              else
              {
                xxTempNewValue += paddingValue * yBasis[b];
                yyTempNewValue += paddingValue * deriv[b];
                zzTempNewValue += paddingValue * yBasis[b];
              }
            } // end b
            grad[0] += xxTempNewValue * zBasis[c];
            grad[1] += yyTempNewValue * zBasis[c];
            grad[2] += zzTempNewValue * deriv[c];
          } // end Z in range
          else
          {
            grad[0] += paddingValue * zBasis[c];
            grad[1] += paddingValue * zBasis[c];
            grad[2] += paddingValue * deriv[c];
          }
        } // end c
      } // end padding value is different from NaN
      else if (previous[0] >= 0.f && previous[0] < (sourceImage->nx - 1) &&
        previous[1] >= 0.f && previous[1] < (sourceImage->ny - 1) &&
        previous[2] >= 0.f && previous[2] < (sourceImage->nz - 1))
      {
        for (c = 0; c < 2; c++)
        {
          Z = previous[2] + c;
          zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
          xxTempNewValue = 0.0;
          yyTempNewValue = 0.0;
          zzTempNewValue = 0.0;
          for (b = 0; b < 2; b++)
          {
            Y = previous[1] + b;
            xyzPointer = &zPointer[Y*sourceImage->nx + previous[0]];
            xTempNewValue = 0.0;
            yTempNewValue = 0.0;
            for (a = 0; a < 2; a++)
            {
              X = previous[0] + a; // \todo This statement should not have any effect optimiser will probably remove this statement anyway
              coeff = *xyzPointer;
              xTempNewValue += coeff * deriv[a];
              yTempNewValue += coeff * xBasis[a];
              xyzPointer++;
            } // end a
            xxTempNewValue += xTempNewValue * yBasis[b];
            yyTempNewValue += yTempNewValue * deriv[b];
            zzTempNewValue += yTempNewValue * yBasis[b];
          } // end b
          grad[0] += xxTempNewValue * zBasis[c];
          grad[1] += yyTempNewValue * zBasis[c];
          grad[2] += zzTempNewValue * deriv[c];
        } // end c
      } // end padding value is NaN
      else grad[0] = grad[1] = grad[2] = 0;

      //check for NaNs in gradient coming from NaNs in floating image
      if (grad[0] != grad[0])
        grad[0] = 0.0;
      if (grad[1] != grad[1])
        grad[1] = 0.0;
      if (grad[2] != grad[2])
        grad[2] = 0.0;

      resultGradientPtrX[index] = (PrecisionType)grad[0];
      resultGradientPtrY[index] = (PrecisionType)grad[1];
      resultGradientPtrZ[index] = (PrecisionType)grad[2];
    }
  }
}




//---------------------------------------
// Transformation::BilinearImageGradient
//---------------------------------------
void Transformation::BilinearImageGradient( nifti_image* sourceImage, nifti_image* resultGradientImage )
{
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)resultGradientImage->nx * resultGradientImage->ny;
  long sourceVoxelNumber = (long)sourceImage->nx * sourceImage->ny;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)resultGradientImage->nx * resultGradientImage->ny;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
#endif

  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultGradientImagePtr = static_cast<PrecisionType *>(resultGradientImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  
  PrecisionType paddingValue = this->warpedPaddingValue;

  mat44 sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = sourceImage->sto_ijk;
  else sourceIJKMatrix = sourceImage->qto_ijk;

  // Iteration over the different volume along the 4th axis
  for (int t = 0; t < resultGradientImage->nt; t++)
  {
    PrecisionType *resultGradientPtrX = &resultGradientImagePtr[2 * t*targetVoxelNumber];
    PrecisionType *resultGradientPtrY = &resultGradientPtrX[targetVoxelNumber];

    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType position[3], xBasis[2], yBasis[2], relative, world[2], grad[2];
    PrecisionType deriv[2];
    deriv[0] = -1;
    deriv[1] = 1;
    PrecisionType coeff, xTempNewValue, yTempNewValue;

    int previous[3], a, b, X, Y;
    PrecisionType *xyPointer;

#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, world, position, previous, xBasis, yBasis, relative, grad, coeff, \
              a, b, X, Y, xyPointer, xTempNewValue, yTempNewValue) \
      shared(sourceIntensity, targetVoxelNumber, sourceVoxelNumber, deriv, \
             deformationFieldPtrX, deformationFieldPtrY, paddingValue, \
             sourceIJKMatrix, sourceImage, resultGradientPtrX, resultGradientPtrY)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {
      grad[0] = 0.0;
      grad[1] = 0.0;

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];

      /* real -> voxel; source space */
      position[0] = world[0] * sourceIJKMatrix.m[0][0] + world[1] * sourceIJKMatrix.m[0][1] +
        sourceIJKMatrix.m[0][3];
      position[1] = world[0] * sourceIJKMatrix.m[1][0] + world[1] * sourceIJKMatrix.m[1][1] +
        sourceIJKMatrix.m[1][3];

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      relative = relative > 0 ? relative : 0;
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      relative = relative > 0 ? relative : 0;
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;

      for (b = 0; b < 2; b++)
      {
        Y = previous[1] + b;
        if (Y > -1 && Y < sourceImage->ny)
        {
          xyPointer = &sourceIntensity[Y*sourceImage->nx + previous[0]];
          xTempNewValue = 0.0;
          yTempNewValue = 0.0;
          for (a = 0; a < 2; a++)
          {
            X = previous[0] + a;
            if (X > -1 && X < sourceImage->nx)
            {
              coeff = *xyPointer;
              xTempNewValue += coeff * deriv[a];
              yTempNewValue += coeff * xBasis[a];
            }
            else
            {
              xTempNewValue += paddingValue * deriv[a];
              yTempNewValue += paddingValue * xBasis[a];
            }
            xyPointer++;
          }
          grad[0] += xTempNewValue * yBasis[b];
          grad[1] += yTempNewValue * deriv[b];
        }
        else
        {
          grad[0] += paddingValue * yBasis[b];
          grad[1] += paddingValue * deriv[b];
        }
      }
      if (grad[0] != grad[0]) grad[0] = 0;
      if (grad[1] != grad[1]) grad[1] = 0;

      resultGradientPtrX[index] = (PrecisionType)grad[0];
      resultGradientPtrY[index] = (PrecisionType)grad[1];
    }
  }
}




//--------------------------------------------
// Transformation::CubicSplineImageGradient3D
//--------------------------------------------
void Transformation::CubicSplineImageGradient3D( nifti_image* sourceImage, nifti_image* resultGradientImage )
{
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)resultGradientImage->nx*resultGradientImage->ny*resultGradientImage->nz;
  long sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)resultGradientImage->nx*resultGradientImage->ny*resultGradientImage->nz;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#endif

  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultGradientImagePtr = static_cast<PrecisionType *>(resultGradientImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[targetVoxelNumber];

  PrecisionType paddingValue = this->warpedPaddingValue;

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  // Iteration over the different volume along the 4th axis
  for (int t = 0; t < resultGradientImage->nt; t++)
  {
    PrecisionType *resultGradientPtrX = &resultGradientImagePtr[3 * t*targetVoxelNumber];
    PrecisionType *resultGradientPtrY = &resultGradientPtrX[targetVoxelNumber];
    PrecisionType *resultGradientPtrZ = &resultGradientPtrY[targetVoxelNumber];

    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    int previous[3], c, Z, b, Y, a;

    PrecisionType xBasis[4], yBasis[4], zBasis[4], xDeriv[4], yDeriv[4], zDeriv[4];
    PrecisionType coeff, position[3], relative, world[3], grad[3];
    PrecisionType xxTempNewValue, yyTempNewValue, zzTempNewValue, xTempNewValue, yTempNewValue;
    PrecisionType *zPointer, *yzPointer, *xyzPointer;
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, world, position, previous, xBasis, yBasis, zBasis, xDeriv, yDeriv, zDeriv, relative, grad, coeff, \
              a, b, c, Y, Z, zPointer, yzPointer, xyzPointer, xTempNewValue, yTempNewValue, xxTempNewValue, yyTempNewValue, zzTempNewValue) \
      shared(sourceIntensity, targetVoxelNumber, sourceVoxelNumber, paddingValue, \
             deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
             sourceIJKMatrix, sourceImage, resultGradientPtrX, resultGradientPtrY, resultGradientPtrZ)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
      grad[2] = 0.0;

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      world[2] = (PrecisionType)deformationFieldPtrZ[index];

      /* real -> voxel; source space */
      reg_mat44_mul( sourceIJKMatrix, world, position );

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      previous[2] = nmm_floorInt( position[2] );

      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      interpolantCubicSpline<PrecisionType>( relative, xBasis, xDeriv );

      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      interpolantCubicSpline<PrecisionType>( relative, yBasis, yDeriv );

      // basis values along the z axis
      relative = position[2] - (PrecisionType)previous[2];
      interpolantCubicSpline<PrecisionType>( relative, zBasis, zDeriv );

      previous[0]--;
      previous[1]--;
      previous[2]--;

      for (c = 0; c < 4; c++)
      {
        Z = previous[2] + c;
        if (-1 < Z && Z < sourceImage->nz)
        {
          zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
          xxTempNewValue = 0.0;
          yyTempNewValue = 0.0;
          zzTempNewValue = 0.0;
          for (b = 0; b < 4; b++)
          {
            Y = previous[1] + b;
            yzPointer = &zPointer[Y*sourceImage->nx];
            if (-1 < Y && Y < sourceImage->ny)
            {
              xyzPointer = &yzPointer[previous[0]];
              xTempNewValue = 0.0;
              yTempNewValue = 0.0;
              for (a = 0; a < 4; a++)
              {
                if (-1 < (previous[0] + a) && (previous[0] + a) < sourceImage->nx)
                {
                  coeff = *xyzPointer;
                  xTempNewValue += coeff * xDeriv[a];
                  yTempNewValue += coeff * xBasis[a];
                } // previous[0]+a in range
                else
                {
                  xTempNewValue += paddingValue * xDeriv[a];
                  yTempNewValue += paddingValue * xBasis[a];
                }
                xyzPointer++;
              } // a
              xxTempNewValue += xTempNewValue * yBasis[b];
              yyTempNewValue += yTempNewValue * yDeriv[b];
              zzTempNewValue += yTempNewValue * yBasis[b];
            } // Y in range
            else
            {
              xxTempNewValue += paddingValue * yBasis[b];
              yyTempNewValue += paddingValue * yDeriv[b];
              zzTempNewValue += paddingValue * yBasis[b];
            }
          } // b
          grad[0] += xxTempNewValue * zBasis[c];
          grad[1] += yyTempNewValue * zBasis[c];
          grad[2] += zzTempNewValue * zDeriv[c];
        } // Z in range
        else
        {
          grad[0] += paddingValue * zBasis[c];
          grad[1] += paddingValue * zBasis[c];
          grad[2] += paddingValue * zDeriv[c];
        }
      } // c

      grad[0] = grad[0] == grad[0] ? grad[0] : 0.0;
      grad[1] = grad[1] == grad[1] ? grad[1] : 0.0;
      grad[2] = grad[2] == grad[2] ? grad[2] : 0.0;
    

      resultGradientPtrX[index] = (PrecisionType)grad[0];
      resultGradientPtrY[index] = (PrecisionType)grad[1];
      resultGradientPtrZ[index] = (PrecisionType)grad[2];
    }
  }
}




//--------------------------------------------
// Transformation::CubicSplineImageGradient2D
//--------------------------------------------
void Transformation::CubicSplineImageGradient2D( nifti_image* sourceImage, nifti_image* resultGradientImage )
{
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)resultGradientImage->nx*resultGradientImage->ny;
  long sourceVoxelNumber = (long)sourceImage->nx * sourceImage->ny;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)resultGradientImage->nx*resultGradientImage->ny;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny;
#endif

  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultGradientImagePtr = static_cast<PrecisionType *>(resultGradientImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  
  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  // Iteration over the different volume along the 4th axis
  for (int t = 0; t < resultGradientImage->nt; t++)
  {
    PrecisionType *resultGradientPtrX = &resultGradientImagePtr[t * 2 * targetVoxelNumber];
    PrecisionType *resultGradientPtrY = &resultGradientPtrX[targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    int previous[2], b, Y, a;
    bool bg;
    PrecisionType xBasis[4], yBasis[4], xDeriv[4], yDeriv[4];
    PrecisionType coeff, position[3], relative, world[3], grad[2];
    PrecisionType xTempNewValue, yTempNewValue;
    PrecisionType *yPointer, *xyPointer;
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, world, position, previous, xBasis, yBasis, xDeriv, yDeriv, relative, grad, coeff, bg, \
              a, b, Y, yPointer, xyPointer, xTempNewValue, yTempNewValue) \
      shared(sourceIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, \
             sourceIJKMatrix, sourceImage, resultGradientPtrX, resultGradientPtrY)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {

      grad[0] = 0.0;
      grad[1] = 0.0;

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];

      /* real -> voxel; source space */
      position[0] = world[0] * sourceIJKMatrix->m[0][0] + world[1] * sourceIJKMatrix->m[0][1] +
        sourceIJKMatrix->m[0][3];
      position[1] = world[0] * sourceIJKMatrix->m[1][0] + world[1] * sourceIJKMatrix->m[1][1] +
        sourceIJKMatrix->m[1][3];

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis, xDeriv );
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis, yDeriv );

      previous[0]--;
      previous[1]--;

      bg = false;
      for (b = 0; b < 4; b++)
      {
        Y = previous[1] + b;
        yPointer = &sourceIntensity[Y*sourceImage->nx];
        if (-1 < Y && Y < sourceImage->ny)
        {
          xyPointer = &yPointer[previous[0]];
          xTempNewValue = 0.0;
          yTempNewValue = 0.0;
          for (a = 0; a < 4; a++)
          {
            if (-1 < (previous[0] + a) && (previous[0] + a) < sourceImage->nx)
            {
              coeff = (PrecisionType)*xyPointer;
              xTempNewValue += coeff * xDeriv[a];
              yTempNewValue += coeff * xBasis[a];
            }
            else bg = true;
            xyPointer++;
          }
          grad[0] += (xTempNewValue * yBasis[b]);
          grad[1] += (yTempNewValue * yDeriv[b]);
        }
        else bg = true;
      }

      if (bg == true)
      {
        grad[0] = 0.0;
        grad[1] = 0.0;
      }
    
      resultGradientPtrX[index] = (PrecisionType)grad[0];
      resultGradientPtrY[index] = (PrecisionType)grad[1];
    }
  }
}




//---------------------------------------------
// Transformation::CheckDVFImageUpdateRequired
//---------------------------------------------
bool Transformation::CheckDVFImageUpdateRequired( nifti_image* targetImageIn )
{
  // No need to check any further if this was already set before
  if ( true == this->dvfImageUpdateRequired )
  {
    return true;
  }

  // Definitively need to calcualte a DVF if it does not exist
  if ( nullptr == this->deformationVectorFieldImage )
  {
    return true;
  }

  if (this->CheckImageGeometryEquality( targetImageIn, this->deformationVectorFieldImage ))
  {
    // If the images have the same basic geometry, no update is required
    return false;
  }
  else
  {
    // Geometries differ, so DVF update is required. 
    return true;
  }
}




//--------------------------------------------
// Transformation::CheckImageGeometryEquality
//--------------------------------------------
bool Transformation::CheckImageGeometryEquality( nifti_image * img1, nifti_image * img2 )
{
  // Check the geometry of img1 and img2
  if (img1->nx != img2->nx) return false;
  if (img1->ny != img2->ny) return false;
  if (img1->nz != img2->nz) return false;
  if (img1->sform_code != img2->sform_code) return false;
  if (img1->qform_code != img2->qform_code) return false;

  // Check sform-matrix or qform matrix otherwise
  if (img1->sform_code > 0)
  {
    // Assuming consistency between sto_ijk and sto_xyz
    if (img1->sto_ijk.m[0][0] != img2->sto_ijk.m[0][0]) return false;
    if (img1->sto_ijk.m[0][1] != img2->sto_ijk.m[0][1]) return false;
    if (img1->sto_ijk.m[0][2] != img2->sto_ijk.m[0][2]) return false;
    if (img1->sto_ijk.m[0][3] != img2->sto_ijk.m[0][3]) return false;

    if (img1->sto_ijk.m[1][0] != img2->sto_ijk.m[1][0]) return false;
    if (img1->sto_ijk.m[1][1] != img2->sto_ijk.m[1][1]) return false;
    if (img1->sto_ijk.m[1][2] != img2->sto_ijk.m[1][2]) return false;
    if (img1->sto_ijk.m[1][3] != img2->sto_ijk.m[1][3]) return false;

    if (img1->sto_ijk.m[2][0] != img2->sto_ijk.m[2][0]) return false;
    if (img1->sto_ijk.m[2][1] != img2->sto_ijk.m[2][1]) return false;
    if (img1->sto_ijk.m[2][2] != img2->sto_ijk.m[2][2]) return false;
    if (img1->sto_ijk.m[2][3] != img2->sto_ijk.m[2][3]) return false;
  }
  else
  {
    // Assuming consistency between qto_ijk and qto_xyz
    if (img1->qto_ijk.m[0][0] != img2->qto_ijk.m[0][0]) return false;
    if (img1->qto_ijk.m[0][1] != img2->qto_ijk.m[0][1]) return false;
    if (img1->qto_ijk.m[0][2] != img2->qto_ijk.m[0][2]) return false;
    if (img1->qto_ijk.m[0][3] != img2->qto_ijk.m[0][3]) return false;

    if (img1->qto_ijk.m[1][0] != img2->qto_ijk.m[1][0]) return false;
    if (img1->qto_ijk.m[1][1] != img2->qto_ijk.m[1][1]) return false;
    if (img1->qto_ijk.m[1][2] != img2->qto_ijk.m[1][2]) return false;
    if (img1->qto_ijk.m[1][3] != img2->qto_ijk.m[1][3]) return false;

    if (img1->qto_ijk.m[2][0] != img2->qto_ijk.m[2][0]) return false;
    if (img1->qto_ijk.m[2][1] != img2->qto_ijk.m[2][1]) return false;
    if (img1->qto_ijk.m[2][2] != img2->qto_ijk.m[2][2]) return false;
    if (img1->qto_ijk.m[2][3] != img2->qto_ijk.m[2][3]) return false;
  }
  
  // Only return true if all of the above were correct
  return true;
}




//---------------------------------------------
// Transformation::CubicSplineTransformImage3D
//---------------------------------------------
void Transformation::ClearDeformationVectorFieldImage()
{
  // Clear DVF image if it exists and set the pointer to null
  if (nullptr != this->deformationVectorFieldImage)
  {
    nifti_image_free( this->deformationVectorFieldImage );
    this->deformationVectorFieldImage = nullptr;;
  }
}



//---------------------------------------------
// Transformation::CubicSplineTransformImage3D
//---------------------------------------------
void Transformation::CubicSplineTransformImage3D( nifti_image* sourceImage, nifti_image* warpedImage )
{
#ifdef _WIN32
  long  index;
  long resultVoxelNumber = (long)warpedImage->nx * warpedImage->ny * warpedImage->nz;
  long sourceVoxelNumber = (long)sourceImage->nx * sourceImage->ny * sourceImage->nz;
#else
  size_t  index;
  size_t resultVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny*warpedImage->nz;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType*>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[resultVoxelNumber];
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[resultVoxelNumber];

  // prepare for openMP loop
  PrecisionType padVal = this->warpedPaddingValue;

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  // Iteration over the different volume along the 4th axis
  for (size_t t = 0; t < (size_t)warpedImage->nt * warpedImage->nu; t++)
  {
    PrecisionType* resultIntensity = &resultIntensityPtr[t*resultVoxelNumber];
    PrecisionType* sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType xBasis[4], yBasis[4], zBasis[4], relative;
    int a, b, c, Y, Z, previous[3];

    PrecisionType *zPointer, *yzPointer, *xyzPointer;
    PrecisionType xTempNewValue, yTempNewValue, intensity, world[3], position[3];
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous, xBasis, yBasis, zBasis, relative, \
              a, b, c, Y, Z, zPointer, yzPointer, xyzPointer, xTempNewValue, yTempNewValue) \
      shared(sourceIntensity, resultIntensity, resultVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
             sourceIJKMatrix, sourceImage, padVal)
#endif // _OPENMP
    for (index = 0; index < resultVoxelNumber; index++)
    {

      intensity = (PrecisionType)(0.0);

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      world[2] = (PrecisionType)deformationFieldPtrZ[index];

      /* real -> voxel; source space */
      reg_mat44_mul( sourceIJKMatrix, world, position );

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      previous[2] = nmm_floorInt( position[2] );

      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );
      // basis values along the z axis
      relative = position[2] - (PrecisionType)previous[2];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, zBasis );

      --previous[0];
      --previous[1];
      --previous[2];

      for (c = 0; c < 4; c++)
      {
        Z = previous[2] + c;
        zPointer = &sourceIntensity[Z * sourceImage->nx * sourceImage->ny];
        yTempNewValue = 0.0;
        for (b = 0; b < 4; b++)
        {
          Y = previous[1] + b;
          yzPointer = &zPointer[Y * sourceImage->nx];
          xyzPointer = &yzPointer[previous[0]];
          xTempNewValue = 0.0;
          for (a = 0; a < 4; a++)
          {
            if (-1 < (previous[0] + a) && (previous[0] + a) < sourceImage->nx && 
                -1 < Z && Z < sourceImage->nz &&
                -1 < Y && Y < sourceImage->ny)
            {
              xTempNewValue += (PrecisionType)*xyzPointer * xBasis[a];
            }
            else
            {
              // paddingValue
              xTempNewValue += padVal * xBasis[a];
            }
            xyzPointer++;
          }
          yTempNewValue += xTempNewValue * yBasis[b];
        }
        intensity += yTempNewValue * zBasis[c];
      }

       resultIntensity[index] = (PrecisionType)intensity;
    }
  }
}




//---------------------------------------------
// Transformation::CubicSplineTransformImage2D
//---------------------------------------------
void Transformation::CubicSplineTransformImage2D( nifti_image* sourceImage, nifti_image* warpedImage ) 
{
  // The resampling scheme is applied along each time
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)warpedImage->nx*warpedImage->ny;
  long sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];

  // prepare for openMP loop
  PrecisionType padVal = this->warpedPaddingValue;


  mat44 sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = sourceImage->sto_ijk;
  else sourceIJKMatrix = sourceImage->qto_ijk;



  for (int t = 0; t < warpedImage->nt*warpedImage->nu; t++)
  {
#ifndef NDEBUG
    printf( "[NiftyReg DEBUG] 2D Cubic spline resampling of volume number %i\n", t );
#endif

    PrecisionType *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType xBasis[4], yBasis[4], relative;
    int a, b, Y, previous[2];

    PrecisionType *yPointer, *xyPointer;
    PrecisionType xTempNewValue, intensity, world[2], position[2];
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous, xBasis, yBasis, relative, \
              a, b, Y, yPointer, xyPointer, xTempNewValue) \
      shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, \
             sourceIJKMatrix, sourceImage, padVal )
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {

      intensity = 0.0;

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      /* real -> voxel; source space */
      position[0] = world[0] * sourceIJKMatrix.m[0][0] + world[1] * sourceIJKMatrix.m[0][1] +
        sourceIJKMatrix.m[0][3];
      position[1] = world[0] * sourceIJKMatrix.m[1][0] + world[1] * sourceIJKMatrix.m[1][1] +
        sourceIJKMatrix.m[1][3];

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );

      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );

      previous[0]--;
      previous[1]--;

      for (b = 0; b < 4; b++)
      {
        Y = previous[1] + b;
        yPointer = &sourceIntensity[Y*sourceImage->nx];
        xyPointer = &yPointer[previous[0]];
        xTempNewValue = 0.0;
        for (a = 0; a < 4; a++)
        {
          if (-1 < (previous[0] + a) && (previous[0] + a) < sourceImage->nx &&
            -1 < Y && Y < sourceImage->ny)
          {
            xTempNewValue += (PrecisionType)*xyPointer * xBasis[a];
          }
          else
          {
            // paddingValue x
            xTempNewValue += padVal * xBasis[a];
          }
          xyPointer++;
        }
        intensity += xTempNewValue * yBasis[b];
      }

      resultIntensity[index] = (PrecisionType)intensity;
    }
  }
}




//--------------------------------------------------
// Transformation::NearestNeighbourTransformImage3D
//--------------------------------------------------
void Transformation::NearestNeighbourTransformImage3D( nifti_image* sourceImage, nifti_image* warpedImage )
{
  // The resampling scheme is applied along each time
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)warpedImage->nx*warpedImage->ny*warpedImage->nz;
  long sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny*warpedImage->nz;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[targetVoxelNumber];

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  for (int t = 0; t < warpedImage->nt*warpedImage->nu; t++)
  {
    PrecisionType *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    // prepare for openMP loop
    PrecisionType padVal = this->warpedPaddingValue;

    PrecisionType intensity;
    PrecisionType world[3];
    PrecisionType position[3];
    int previous[3];

#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous) \
      shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
             sourceIJKMatrix, sourceImage, padVal)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {
      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      world[2] = (PrecisionType)deformationFieldPtrZ[index];

      /* real -> voxel; source space */
      reg_mat44_mul( sourceIJKMatrix, world, position );

      previous[0] = (int)std::round( position[0] );
      previous[1] = (int)std::round( position[1] );
      previous[2] = (int)std::round( position[2] );

      if (-1 < previous[2] && previous[2] < sourceImage->nz &&
          -1 < previous[1] && previous[1] < sourceImage->ny &&
          -1 < previous[0] && previous[0] < sourceImage->nx)
      {
        intensity = sourceIntensity[(previous[2] * sourceImage->ny + previous[1]) * sourceImage->nx + previous[0]];
        resultIntensity[index] = intensity;
      }
      else resultIntensity[index] = padVal;
    }
  }
}



//--------------------------------------------------
// Transformation::NearestNeighbourTransformImage2D
//--------------------------------------------------
void Transformation::NearestNeighbourTransformImage2D( nifti_image* sourceImage, nifti_image* warpedImage )
{
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)warpedImage->nx*warpedImage->ny;
  long sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny;
#else
  size_t  index;
  size_t targetVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];

  // prepare for openMP loop
  PrecisionType padVal = this->warpedPaddingValue;
  
  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  for (int t = 0; t < warpedImage->nt * warpedImage->nu; t++)
  {
    PrecisionType *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType intensity;
    PrecisionType world[2];
    PrecisionType position[2];
    int previous[2];

#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous) \
      shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, \
             sourceIJKMatrix, sourceImage, padVal)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {

      world[0] = (PrecisionType)deformationFieldPtrX[index];
      world[1] = (PrecisionType)deformationFieldPtrY[index];
      /* real -> voxel; source space */
      position[0] = world[0] * sourceIJKMatrix->m[0][0] + world[1] * sourceIJKMatrix->m[0][1] + sourceIJKMatrix->m[0][3];
      position[1] = world[0] * sourceIJKMatrix->m[1][0] + world[1] * sourceIJKMatrix->m[1][1] + sourceIJKMatrix->m[1][3];

      previous[0] = static_cast<int>(std::round( position[0] ));
      previous[1] = static_cast<int>(std::round( position[1] ));

      if (-1 < previous[1] && previous[1] < sourceImage->ny &&
          -1 < previous[0] && previous[0] < sourceImage->nx)
      {
        intensity = sourceIntensity[previous[1] * sourceImage->nx + previous[0]];
        resultIntensity[index] = intensity;
      }
      else resultIntensity[index] = (PrecisionType) padVal;
    }
  }
}



//-----------------------------------------
// Transformation::TrilinearTransformImage
//-----------------------------------------
void Transformation::TrilinearTransformImage( nifti_image* sourceImage, nifti_image* warpedImage )
{
  // The resampling scheme is applied along each time
#ifdef _WIN32
  long index;
  long targetVoxelNumber = (long)warpedImage->nx*warpedImage->ny*warpedImage->nz;
  long sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#else
  size_t index;
  size_t targetVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny*warpedImage->nz;
  size_t sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny*sourceImage->nz;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[targetVoxelNumber];

  // prepare for openMP loop
  PrecisionType padVal = this->warpedPaddingValue;

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  for (int t = 0; t < warpedImage->nt*warpedImage->nu; t++)
  {
    PrecisionType *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType xBasis[2], yBasis[2], zBasis[2], relative;
    int a, b, c, X, Y, Z, previous[3];

    PrecisionType *zPointer, *xyzPointer;
    PrecisionType xTempNewValue, yTempNewValue, intensity, world[3], position[3];
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous, xBasis, yBasis, zBasis, relative, \
              a, b, c, X, Y, Z, zPointer, xyzPointer, xTempNewValue, yTempNewValue) \
      shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
             sourceIJKMatrix, sourceImage, padVal)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; index++)
    {

      intensity = padVal;

      //check in mask and def field not NaN
      if ( deformationFieldPtrX[index] == deformationFieldPtrX[index]
        && deformationFieldPtrY[index] == deformationFieldPtrY[index]
        && deformationFieldPtrZ[index] == deformationFieldPtrZ[index])
      {

        intensity = 0;

        world[0] = (PrecisionType)deformationFieldPtrX[index];
        world[1] = (PrecisionType)deformationFieldPtrY[index];
        world[2] = (PrecisionType)deformationFieldPtrZ[index];

        /* real -> voxel; source space */
        reg_mat44_mul( sourceIJKMatrix, world, position );

        previous[0] = nmm_floorInt( position[0] );
        previous[1] = nmm_floorInt( position[1] );
        previous[2] = nmm_floorInt( position[2] );

        // basis values along the x axis
        relative = position[0] - (PrecisionType)previous[0];
        xBasis[0] = (PrecisionType)(1.0 - relative);
        xBasis[1] = relative;
        // basis values along the y axis
        relative = position[1] - (PrecisionType)previous[1];
        yBasis[0] = (PrecisionType)(1.0 - relative);
        yBasis[1] = relative;
        // basis values along the z axis
        relative = position[2] - (PrecisionType)previous[2];
        zBasis[0] = (PrecisionType)(1.0 - relative);
        zBasis[1] = relative;

        // For efficiency reason two interpolations are here, with and without using a padding value
        if (padVal == padVal)
        {
          // Interpolation using the padding value
          for (c = 0; c < 2; c++)
          {
            Z = previous[2] + c;
            if (Z > -1 && Z < sourceImage->nz)
            {
              zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
              yTempNewValue = 0.0;
              for (b = 0; b < 2; b++)
              {
                Y = previous[1] + b;
                if (Y > -1 && Y < sourceImage->ny)
                {
                  xyzPointer = &zPointer[Y*sourceImage->nx + previous[0]];
                  xTempNewValue = 0.0;
                  for (a = 0; a < 2; a++)
                  {
                    X = previous[0] + a;
                    if (X > -1 && X < sourceImage->nx)
                    {
                      xTempNewValue += *xyzPointer * xBasis[a];
                    } // X
                    else xTempNewValue += padVal * xBasis[a];
                    xyzPointer++;
                  } // a
                  yTempNewValue += xTempNewValue * yBasis[b];
                } // Y
                else yTempNewValue += padVal * yBasis[b];
              } // b
              intensity += yTempNewValue * zBasis[c];
            } // Z
            else intensity += padVal * zBasis[c];
          } // c
        } // padding value is defined
        else if (previous[0] >= 0.f && previous[0] < (sourceImage->nx - 1) &&
          previous[1] >= 0.f && previous[1] < (sourceImage->ny - 1) &&
          previous[2] >= 0.f && previous[2] < (sourceImage->nz - 1))
        {
          for (c = 0; c < 2; c++)
          {
            Z = previous[2] + c;
            zPointer = &sourceIntensity[Z*sourceImage->nx*sourceImage->ny];
            yTempNewValue = 0.0;
            for (b = 0; b < 2; b++)
            {
              Y = previous[1] + b;
              xyzPointer = &zPointer[Y*sourceImage->nx + previous[0]];
              xTempNewValue = 0.0;
              for (a = 0; a < 2; a++)
              {
                X = previous[0] + a;
                xTempNewValue += *xyzPointer * xBasis[a];
                xyzPointer++;
              } // a
              yTempNewValue += xTempNewValue * yBasis[b];
            } // b
            intensity += yTempNewValue * zBasis[c];
          } // c
        } // padding value is not defined
        // The voxel is outside of the source space and thus set to NaN here
        else
        {
          intensity = padVal;
        }
      } // voxel is in the mask

      resultIntensity[index] = (PrecisionType)intensity;
    }
  }
}

//----------------------------------------
// Transformation::BilinearTransformImage
//----------------------------------------
void Transformation::BilinearTransformImage( nifti_image* sourceImage, nifti_image* warpedImage )
{
  // The resampling scheme is applied along each time
#ifdef _WIN32
  long  index;
  long  targetVoxelNumber = (long)warpedImage->nx*warpedImage->ny;
  long  sourceVoxelNumber = (long)sourceImage->nx*sourceImage->ny;
#else
  size_t  index;
  size_t  targetVoxelNumber = (size_t)warpedImage->nx*warpedImage->ny;
  size_t  sourceVoxelNumber = (size_t)sourceImage->nx*sourceImage->ny;
#endif
  PrecisionType *sourceIntensityPtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *resultIntensityPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[targetVoxelNumber];

  // prepare for openMP loop
  PrecisionType padVal = this->warpedPaddingValue;

  mat44 *sourceIJKMatrix;
  if (sourceImage->sform_code > 0)
    sourceIJKMatrix = &(sourceImage->sto_ijk);
  else sourceIJKMatrix = &(sourceImage->qto_ijk);

  for (int t = 0; t < warpedImage->nt*warpedImage->nu; t++)
  {

#ifndef NDEBUG
    printf( "[NiftyReg DEBUG] 2D linear resampling of volume number %i\n", t );
#endif

    PrecisionType *resultIntensity = &resultIntensityPtr[t*targetVoxelNumber];
    PrecisionType *sourceIntensity = &sourceIntensityPtr[t*sourceVoxelNumber];

    PrecisionType xBasis[2], yBasis[2], relative;
    int a, b, X, Y, previous[3];

    PrecisionType *xyPointer;
    PrecisionType xTempNewValue, intensity, world[2], position[2];
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
      private(index, intensity, world, position, previous, xBasis, yBasis, relative, \
              a, b, X, Y, xyPointer, xTempNewValue) \
      shared(sourceIntensity, resultIntensity, targetVoxelNumber, sourceVoxelNumber, \
             deformationFieldPtrX, deformationFieldPtrY, \
             sourceIJKMatrix, sourceImage, padVal)
#endif // _OPENMP
    for (index = 0; index < targetVoxelNumber; ++index)
    {
      intensity = 0;

      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];

      /* real -> voxel; source space */
      position[0] = world[0] * sourceIJKMatrix->m[0][0] + world[1] * sourceIJKMatrix->m[0][1] + sourceIJKMatrix->m[0][3];
      position[1] = world[0] * sourceIJKMatrix->m[1][0] + world[1] * sourceIJKMatrix->m[1][1] + sourceIJKMatrix->m[1][3];

      previous[0] = nmm_floorInt( position[0] );
      previous[1] = nmm_floorInt( position[1] );
      
      // basis values along the x axis
      relative = position[0] - (PrecisionType)previous[0];
      relative = relative > 0 ? relative : 0;
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      
      // basis values along the y axis
      relative = position[1] - (PrecisionType)previous[1];
      relative = relative > 0 ? relative : 0;
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;

      for (b = 0; b < 2; b++)
      {
        Y = previous[1] + b;
        if (Y > -1 && Y < sourceImage->ny)
        {
          xyPointer = &sourceIntensity[Y*sourceImage->nx + previous[0]];
          xTempNewValue = 0.0;
      
          for (a = 0; a < 2; a++)
          {
            X = previous[0] + a;
            if (X > -1 && X < sourceImage->nx)
            {
              xTempNewValue += *xyPointer * xBasis[a];
            }
            else xTempNewValue += padVal * xBasis[a];
            xyPointer++;
          } // a

          intensity += xTempNewValue * yBasis[b];
        } // Y outside
        else intensity += padVal * yBasis[b];
      } // b
      
      resultIntensity[index] = (PrecisionType)intensity;
    }
  }
}




//----------------------------------------------------
// Transformation::CubicSplineTransformImageAdjoint3D
//----------------------------------------------------
void Transformation::CubicSplineTransformImageAdjoint3D( 
  nifti_image * sourceImage, 
  nifti_image * sourceWeightsImage, 
  nifti_image * warpedImage, 
  nifti_image * warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];

  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[4], yBasis[4], zBasis[4], relative, world[3], position[3], *warpedWeiZPtr;
  PrecisionType *warpedWeiXYZPtr, yzBasis, xyzBasis, intensity;
  int index, a, b, c, firstVox[3], Y, Z;
  PrecisionType *warpedImgZPtr, *warpedImgXYZPtr, warpedImgValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, zBasis, relative, world, position, firstVox, \
    a, b, c, Y, Z, warpedImgZPtr, warpedWeiZPtr, warpedWeiXYZPtr, warpedImgXYZPtr, \
    yzBasis, xyzBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      
      //round down to voxel before def field position
      firstVox[0] = nmm_floorInt( position[0] );
      firstVox[1] = nmm_floorInt( position[1] );
      firstVox[2] = nmm_floorInt( position[2] );
      
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );
      
      // basis values along the z axis
      relative = position[2] - (PrecisionType)firstVox[2];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, zBasis );
      
      //as using cubic spline need two voxels before and two after def field position
      //so shift firstVox one more in -ve direction
      --firstVox[0];
      --firstVox[1];
      --firstVox[2];

      //loop over 4 voxels on z axis
      for (c = 0; c < 4; c++)
      {
        Z = firstVox[2] + c;
        if (-1 < Z && Z < warpedImage->nz)
        {
          warpedImgZPtr = &warpedImgPtr[Z * warpedImage->nx * warpedImage->ny];
          warpedWeiZPtr = &warpedWeightsPtr[Z * warpedWeightsImage->nx * warpedWeightsImage->ny];

          //loop over 4 voxels on y axis
          for (b = 0; b < 4; b++)
          {
            Y = firstVox[1] + b;
            if (-1 < Y && Y < warpedImage->ny)
            {
              warpedImgXYZPtr = &warpedImgZPtr[Y * warpedImage->nx + firstVox[0]];
              warpedWeiXYZPtr = &warpedWeiZPtr[Y * warpedWeightsImage->nx + firstVox[0]];
              yzBasis = yBasis[b] * zBasis[c];

              //loop over 4 voxels on x axis
              for (a = 0; a < 4; a++)
              {
                if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
                {
                  xyzBasis = yzBasis * xBasis[a];
                  //#pragma omp atomic
                  warpedWeiXYZPtr[a] += sourceWeightsPtr[index] * xyzBasis;
                  intensity = ((PrecisionType)sourcePtr[index]) * xyzBasis;
                  
                  /// \todo: Decide if original switch statement is necessary at all. 
                  //switch (sourceImage->datatype)
                  //{
                  //case NIFTI_TYPE_FLOAT32:
                  //case NIFTI_TYPE_FLOAT64:
                  warpedImgValue = (PrecisionType)intensity;
                  //  break;
                  //case NIFTI_TYPE_UINT8:
                  //case NIFTI_TYPE_UINT16:
                  //case NIFTI_TYPE_UINT32:
                  //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
                  //  break;
                  //default:
                  //  warpedImgValue = (PrecisionType)reg_round( intensity );
                  //}//switch (dynamicImage->datatype)
//#pragma omp atomic
                  warpedImgXYZPtr[a] += warpedImgValue;
                }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
              }//for (a = 0; a < 4; a++)
            }//if (-1 < Y && Y < warpedDynImage->ny)
          }//for (b = 0; b < 4; b++)
        }//if (-1 < Z && Z < warpedDynImage->nz)
      }//for (c = 0; c < 4; c++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//----------------------------------------------------
// Transformation::CubicSplineTransformImageAdjoint2D
//----------------------------------------------------
void Transformation::CubicSplineTransformImageAdjoint2D( 
  nifti_image * sourceImage, 
  nifti_image * sourceWeightsImage, 
  nifti_image * warpedImage, 
  nifti_image * warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image
 
  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[4], yBasis[4], relative, world[2], position[2], *warpedWeiXYPtr;
  PrecisionType xyBasis, intensity;
  int index, a, b, firstVox[2], Y;
  PrecisionType *warpedImgXYPtr, warpedImgValue;

/*  #if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, relative, world, position, firstVox, \
    a, b, Y, warpedWeiXYPtr, warpedImgXYPtr, xyBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, warpedIJKMatrix, \
    warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP */

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] + (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] + (*warpedIJKMatrix).m[1][3];
      
      //round down to voxel before def field position
      firstVox[0] = nmm_floorInt( position[0] );
      firstVox[1] = nmm_floorInt( position[1] );
      
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );
      
      //as using cubic spline need two voxels before and two after def field position
      //so shift firstVox one more in -ve direction
      --firstVox[0];
      --firstVox[1];

      //loop over 4 voxels on y axis
      for (b = 0; b < 4; b++)
      {
        Y = firstVox[1] + b;
        if (-1 < Y && Y < warpedImage->ny)
        {
          warpedImgXYPtr = &warpedImgPtr[Y * warpedImage->nx + firstVox[0]];
          warpedWeiXYPtr = &warpedWeightsPtr[Y * warpedWeightsImage->nx + firstVox[0]];

          //loop over 4 voxels on x axis
          for (a = 0; a < 4; a++)
          {
            if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
            {
              xyBasis = yBasis[b] * xBasis[a];
              //#pragma omp atomic
              warpedWeiXYPtr[a] += sourceWeightsPtr[index] * xyBasis;
              intensity = ((PrecisionType)sourcePtr[index]) * xyBasis;
              /// \todo: Decide if original switch statement is necessary at all. 
              //switch (sourceImage->datatype)
              //{
              //case NIFTI_TYPE_FLOAT32:
              //case NIFTI_TYPE_FLOAT64:
              warpedImgValue = (PrecisionType)intensity;
              //  break;
              //case NIFTI_TYPE_UINT8:
              //case NIFTI_TYPE_UINT16:
              //case NIFTI_TYPE_UINT32:
              //  warpedImgValue = (ImageTYPE)(intensity > 0 ? reg_round( intensity ) : 0);
              //  break;
              //default:
              //  warpedImgValue = (ImageTYPE)reg_round( intensity );
              //}//switch (dynamicImage->datatype)
//#pragma omp atomic
              warpedImgXYPtr[a] += warpedImgValue;
            }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
          }//for (a = 0; a < 4; a++)
        }//if (-1 < Y && Y < warpedDynImage->ny)
      }//for (b = 0; b < 4; b++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//---------------------------------------------------------
// Transformation::NearestNeighbourTransformImageAdjoint3D
//---------------------------------------------------------
void Transformation::NearestNeighbourTransformImageAdjoint3D( 
  nifti_image * sourceImage, 
  nifti_image * sourceWeightsImage, 
  nifti_image * warpedImage, 
  nifti_image * warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of source image not warped image
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];

  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType world[3], position[3];
  int index, nearestVox[3], warpedInd;

/*  #if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, world, position, nearestVox, warpedInd) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP */

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];
      //now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      //round to nearest voxel to def field position
      nearestVox[0] = static_cast<int>(reg_round( position[0] ));
      nearestVox[1] = static_cast<int>(reg_round( position[1] ));
      nearestVox[2] = static_cast<int>(reg_round( position[2] ));

      //if nearestVox in image push value and weight
      if (-1 < nearestVox[0] && nearestVox[0] < warpedImage->nx &&
          -1 < nearestVox[1] && nearestVox[1] < warpedImage->ny &&
          -1 < nearestVox[2] && nearestVox[2] < warpedImage->nz)
      {
        warpedInd = (nearestVox[2] * warpedImage->ny + nearestVox[1]) * warpedImage->nx + nearestVox[0];
        //#pragma omp atomic
        warpedImgPtr[warpedInd] += sourcePtr[index];
        //#pragma omp atomic
        warpedWeightsPtr[warpedInd] += sourceWeightsPtr[index];
      }
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//---------------------------------------------------------
// Transformation::NearestNeighbourTransformImageAdjoint2D
//---------------------------------------------------------
void Transformation::NearestNeighbourTransformImageAdjoint2D( 
  nifti_image * sourceImage, 
  nifti_image * sourceWeightsImage, 
  nifti_image * warpedImage, 
  nifti_image * warpedWeightsImage )
{
  //some useful pointer definitions
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image

  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType world[2], position[2];
  int index, nearestVox[2], warpedInd;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, world, position, nearestVox, warpedInd) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP */

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] + (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] + (*warpedIJKMatrix).m[1][3];
      
      //round to nearest voxel to def field position
      nearestVox[0] = static_cast<int>(reg_round( position[0] ));
      nearestVox[1] = static_cast<int>(reg_round( position[1] ));

      //if nearestVox in image push value and weight
      if (-1 < nearestVox[0] && nearestVox[0] < warpedImage->nx &&
          -1 < nearestVox[1] && nearestVox[1] < warpedImage->ny)
      {
        warpedInd = nearestVox[1] * warpedImage->nx + nearestVox[0];
        //#pragma omp atomic
        warpedImgPtr[warpedInd] += sourcePtr[index];
        //#pragma omp atomic
        warpedWeightsPtr[warpedInd] += sourceWeightsPtr[index];
      }
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//------------------------------------------------
// Transformation::TrilinearTransformImageAdjoint
//------------------------------------------------
void Transformation::TrilinearTransformImageAdjoint( 
  nifti_image* sourceImage, 
  nifti_image* sourceWeightsImage, 
  nifti_image* warpedImage, 
  nifti_image* warpedWeightsImage )
{
  //some useful pointer definitions
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of source image not warped image
  PrecisionType* deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];
  
  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0) 
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[2], yBasis[2], zBasis[2], relative, world[3], position[3], *warpedWeiZPtr;
  PrecisionType *warpedWeiXYZPtr, yzBasis, xyzBasis, intensity;
  int index, a, b, c, firstVox[3], Y, Z;
  PrecisionType *warpedImgZPtr, *warpedImgXYZPtr, warpedImgValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, zBasis, relative, world, position, firstVox, \
    a, b, c, Y, Z, warpedImgZPtr, warpedWeiZPtr, warpedWeiXYZPtr, warpedImgXYZPtr, \
    yzBasis, xyzBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP */

  // as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; ++index)
  {
    // only push value if not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      // find where the def field maps this source voxel to in the warped image
      // first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];

      // now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      
      // round down to voxel before def field position
      firstVox[0] = nmm_floorInt( position[0] );
      firstVox[1] = nmm_floorInt( position[1] );
      firstVox[2] = nmm_floorInt( position[2] );
      
      // get basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;
      // basis values along the z axis
      relative = position[2] - (PrecisionType)firstVox[2];
      zBasis[0] = (PrecisionType)(1.0 - relative);
      zBasis[1] = relative;

      //loop over 2 voxels on z axis
      for (c = 0; c < 2; c++)
      {
        Z = firstVox[2] + c;
        if (-1 < Z && Z < warpedImage->nz)
        {
          warpedImgZPtr = &warpedImgPtr[Z * warpedImage->nx * warpedImage->ny];
          warpedWeiZPtr = &warpedWeightsPtr[Z * warpedWeightsImage->nx * warpedWeightsImage->ny];

          //loop over 2 voxels on y axis
          for (b = 0; b < 2; b++)
          {
            Y = firstVox[1] + b;
            if (-1 < Y && Y < warpedImage->ny)
            {
              warpedImgXYZPtr = &warpedImgZPtr[Y * warpedImage->nx + firstVox[0]];
              warpedWeiXYZPtr = &warpedWeiZPtr[Y * warpedWeightsImage->nx + firstVox[0]];
              yzBasis = yBasis[b] * zBasis[c];

              //loop over 2 voxels on x axis
              for (a = 0; a < 2; a++)
              {
                if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
                {
                  xyzBasis = yzBasis * xBasis[a];
                  
                  //#pragma omp atomic
                  warpedWeiXYZPtr[a] += sourceWeightsPtr[index] * xyzBasis;
                  intensity = ((PrecisionType)sourcePtr[index]) * xyzBasis;
                  /// \todo: Decide if original switch statement is necessary at all. 
                  //switch (sourceImage->datatype)
                  //{
                  //case NIFTI_TYPE_FLOAT32:
                  //case NIFTI_TYPE_FLOAT64:
                  warpedImgValue = (PrecisionType)intensity;
                  //  break;
                  //case NIFTI_TYPE_UINT8:
                  //case NIFTI_TYPE_UINT16:
                  //case NIFTI_TYPE_UINT32:
                  //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
                  //  break;
                  //default:
                  //  warpedImgValue = (PrecisionType)reg_round( intensity );
                  //}//switch (dynamicImage->datatype)
//#pragma omp atomic
                  warpedImgXYZPtr[a] += warpedImgValue;
                }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
              }//for (a = 0; a < 2; a++)
            }//if (-1 < Y && Y < warpedDynImage->ny)
          }//for (b = 0; b < 2; b++)
        }//if (-1 < Z && Z < warpedDynImage->nz)
      }//for (c = 0; c < 2; c++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//-----------------------------------------------
// Transformation::BilinearTransformImageAdjoint
//-----------------------------------------------
void Transformation::BilinearTransformImageAdjoint( 
  nifti_image* sourceImage, 
  nifti_image* sourceWeightsImage, 
  nifti_image* warpedImage, 
  nifti_image* warpedWeightsImage )
{
  //some useful pointers
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image

  //matrix from real space to warped image space
  mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[2], yBasis[2], relative, world[2], position[2], *warpedWeiXYPtr;
  PrecisionType xyBasis, intensity;
  int index, a, b, firstVox[2], Y;
  PrecisionType *warpedImgXYPtr, warpedDynValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, relative, world, position, firstVox, \
    a, b, Y, warpedWeiXYPtr, warpedImgXYPtr, xyBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, maskPtr, warpedIJKMatrix, \
    warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP */

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] +
        (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] +
        (*warpedIJKMatrix).m[1][3];
      //round down to voxel before def field position
      firstVox[0] = nmm_floorInt( position[0] );
      firstVox[1] = nmm_floorInt( position[1] );
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;

      //loop over 2 voxels on y axis
      for (b = 0; b < 2; b++)
      {
        Y = firstVox[1] + b;
        if (-1 < Y && Y < warpedImage->ny)
        {
          warpedImgXYPtr = &warpedImgPtr[Y*warpedImage->nx + firstVox[0]];
          warpedWeiXYPtr = &warpedWeightsPtr[Y*warpedWeightsImage->nx + firstVox[0]];

          //loop over 2 voxels on x axis
          for (a = 0; a < 2; a++)
          {
            if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
            {
              xyBasis = yBasis[b] * xBasis[a];
              //#pragma omp atomic
              warpedWeiXYPtr[a] += sourceWeightsPtr[index] * xyBasis;
              intensity = ((PrecisionType)sourcePtr[index]) * xyBasis;
              /// \todo: Decide if original switch statement is necessary at all. 
              /// \todo: Remove reg_round macro with inline function
              //switch (sourceImage->datatype)
              //{
              //case NIFTI_TYPE_FLOAT32:
              //case NIFTI_TYPE_FLOAT64:
              warpedDynValue = (PrecisionType)intensity;
              //  break;
              //case NIFTI_TYPE_UINT8:
              //case NIFTI_TYPE_UINT16:
              //case NIFTI_TYPE_UINT32:
              //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
              //  break;
              //default:
              //  warpedImgValue = (PrecisionType)reg_round( intensity );
              //}//switch (dynamicImage->datatype)
//#pragma omp atomic
              warpedImgXYPtr[a] += warpedDynValue;
            }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
          }//for (a = 0; a < 2; a++)
        }//if (-1 < Y && Y < warpedDynImage->ny)
      }//for (b = 0; b < 2; b++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}



//----------------------------------------
// Transformation::ReorientateVectorImage
//----------------------------------------
void Transformation::ReorientateVectorImage(nifti_image* vectorFieldImageToReorientate, mat44 reorientationMatrix)
{
	PrecisionType *vecPtrX = static_cast<PrecisionType *>(vectorFieldImageToReorientate->data);
	PrecisionType *vecPtrY = &vecPtrX[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];
	PrecisionType *vecPtrZ = NULL;
	if (vectorFieldImageToReorientate->nu == 3)
		vecPtrZ = &vecPtrY[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];

	for (int n = 0; n < vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt; n++)
	{
		PrecisionType reorientatedValues[3] = { 0.0, 0.0, 0.0 };
		if (vecPtrZ == NULL) // 2D
		{
			reorientatedValues[0] =
				reorientationMatrix.m[0][0] * vecPtrX[n] +
				reorientationMatrix.m[1][0] * vecPtrY[n];
			reorientatedValues[1] =
				reorientationMatrix.m[0][1] * vecPtrX[n] +
				reorientationMatrix.m[1][1] * vecPtrY[n];
			vecPtrX[n] = reorientatedValues[0];
			vecPtrY[n] = reorientatedValues[1];
		}
		else // 3D
		{
			reorientatedValues[0] =
				reorientationMatrix.m[0][0] * vecPtrX[n] +
				reorientationMatrix.m[1][0] * vecPtrY[n] +
				reorientationMatrix.m[2][0] * vecPtrZ[n];
			reorientatedValues[1] =
				reorientationMatrix.m[0][1] * vecPtrX[n] +
				reorientationMatrix.m[1][1] * vecPtrY[n] +
				reorientationMatrix.m[2][1] * vecPtrZ[n];
			reorientatedValues[2] =
				reorientationMatrix.m[0][2] * vecPtrX[n] +
				reorientationMatrix.m[1][2] * vecPtrY[n] +
				reorientationMatrix.m[2][2] * vecPtrZ[n];
			vecPtrX[n] = reorientatedValues[0];
			vecPtrY[n] = reorientatedValues[1];
			vecPtrZ[n] = reorientatedValues[2];
		}
	}

	return;
}

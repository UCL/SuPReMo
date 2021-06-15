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

#include "Supremo.h"


/** Test 4x4 matrices for equality element by element. Matrix size of 4x4 is assumed. 
 * \param mat1 First matrix
 * \param mat2 Second matrix
 */
int checkMat44Equality( const mat44& mat1, const mat44& mat2 )
{
  for (unsigned int i = 0; i < 4; ++i)
  {
    for (unsigned int j = 0; j < 4; ++j)
    {
      if (mat1.m[i][j] != mat2.m[i][j])
        return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}



/** Test images for equal geometry
 *  The test will complete all test before returning EXIT_SUCCESS or EXIT_FALIURE
 *  Items tested are:
 *   - nx, ny, nz
 *   - dx, dy, dz
 *   - sform_code, qform_code
 *   - sto_ijk, qto_ijk, sto_xyz, qto_xyz
 * sform and qform matrices are only checked if sform and qform codes are unequal to zero respectively. 
 *
 * \param imageToTest Pointer to the nifti image that is compared to the expected outcome
 * \param expectedImage The expected outcome image
 */
int checkImageGeometryEquality( nifti_image* imageToTest, nifti_image* expectedImage )
{
   bool imagesEqual = true;

  // Check basic image parameters 
  // size (nx/ny/xz)
  if (imageToTest->nx != expectedImage->nx)
  {
    supremo_print_error( "Generated image and requested image differ in x-dimension" );
    imagesEqual = false;
  }
  if (imageToTest->ny != expectedImage->ny)
  {
    supremo_print_error( "Generated image and requested image differ in y-dimension" );
    imagesEqual = false;
  }
  if (imageToTest->nz != expectedImage->nz)
  {
    supremo_print_error( "Generated image and requested image differ in z-dimension" );
    imagesEqual = false;
  }
  // spacing (dx/dy/dz)
  if (imageToTest->dx != expectedImage->dx)
  {
    supremo_print_error( "Generated image and requested image differ in x-spacing" );
    imagesEqual = false;
  }
  if (imageToTest->dy != expectedImage->dy)
  {
    supremo_print_error( "Generated image and requested image differ in y-spacing" );
    imagesEqual = false;
  }
  if (imageToTest->dz != expectedImage->dz)
  {
    supremo_print_error( "Generated image and requested image differ in z-spacing" );
    imagesEqual = false;
  }

  // Check geometry i.e. sform and qform 
  if (imageToTest->sform_code != expectedImage->sform_code)
  {
    supremo_print_error( "Generated image and requested image differ in sform_code" );
    imagesEqual = false;
  }
  if (imageToTest->qform_code != expectedImage->qform_code)
  {
    supremo_print_error( "Generated image and requested image differ in sform_code" );
    imagesEqual = false;
  }
  if (imageToTest->sform_code != 0)
  {
    if (checkMat44Equality( imageToTest->sto_ijk, expectedImage->sto_ijk ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated image and requested image differ in sto_ijk" );
      imagesEqual = false;
    }
    if (checkMat44Equality( imageToTest->sto_xyz, expectedImage->sto_xyz ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated image and requested image differ in sto_xyz" );
      imagesEqual = false;
    }
  }
  if (imageToTest->qform_code != 0)
  {
    if (checkMat44Equality( imageToTest->qto_ijk, expectedImage->qto_ijk ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated image and requested image differ in qto_ijk" );
      imagesEqual = false;
    }
    if (checkMat44Equality( imageToTest->qto_xyz, expectedImage->qto_xyz ) == EXIT_FAILURE)
    {
      supremo_print_error( "Generated image and requested image differ in qto_xyz" );
      imagesEqual = false;
    }
  }

  // Return the correct return value
  if (imagesEqual)
  {
    return EXIT_SUCCESS;
  }
  else return EXIT_FAILURE;
}




/** Test images for equality
 *  The test will complete the image geometry test and in addition will test the contents. Currently only 
 *  three dimensional images are supported. 
 *  Items tested are:
 *   - nbyper
 *   - voxel-by-voxel comparision
 * \param imageToTest Pointer to the nifti image that is compared to the expected outcome.
 * \param expectedImage The expected outcome image. 
 * \param allowedVoxelDeviation Allow for numerical inaccuracies. 
 * \param ignoreBoundaySize Ignore differences at the border for the given size. 
 * \param maxAllowedNumberOfDeviations Allow a certain deviations up to the given number. 
 */
template<class VoxelType>
int checkImageEquality( nifti_image* imageToTest, 
  nifti_image* expectedImage, 
  VoxelType allowedVoxelDeviation, 
  unsigned int ignoreBoundaySize = 0, 
  unsigned int maxAllowedNumberOfDeviations = 0)
{
  // Track the equality results
  bool imagesEqual = true;
 
  // Check the image geometry first
  if (checkImageGeometryEquality( imageToTest, expectedImage ) == EXIT_FAILURE) imagesEqual = false;
  
  // Check the type
  if (imageToTest->nbyper != expectedImage->nbyper)
  {
    supremo_print_error( "Generated image and requested image differ in nbyper" );
    imagesEqual = false;
  }

  if (sizeof( VoxelType ) != expectedImage->nbyper)
  {
    supremo_print_error( "Wrong voxel type assumed" );
    imagesEqual = false;
  }

  // Check image contents
  VoxelType* imageDataPter = (VoxelType*)(imageToTest)->data;
  VoxelType* expectedImageDataPter = (VoxelType*)(expectedImage)->data;
  unsigned int numberOfDeviations = 0;
  VoxelType maxDeviation = 0;
  if (ignoreBoundaySize != 0)
  {
    // Make sure below is only executed when 3D
    if (expectedImage->dim[0] == 3)
    {
      for (unsigned int z = ignoreBoundaySize; z < expectedImage->nz - ignoreBoundaySize; ++z)
      {
        for (unsigned int y = ignoreBoundaySize; y < expectedImage->ny - ignoreBoundaySize; ++y)
        {
          for (unsigned int x = ignoreBoundaySize; x < expectedImage->nz - ignoreBoundaySize; ++x)
          {
            int idx = (expectedImage->ny * z + y) * expectedImage->nx + x;
            VoxelType absDif = std::abs( imageDataPter[idx] - expectedImageDataPter[idx] );

            if (absDif > allowedVoxelDeviation)
            {
              // Make sure that if both images are not a number, then do something else
              if (imageDataPter[idx] != imageDataPter[idx] || expectedImageDataPter[idx] != expectedImageDataPter[idx]) continue;

              maxDeviation = maxDeviation < absDif ? absDif : maxDeviation;
              ++numberOfDeviations;
            }
          }
        }
      }
    }
    if (expectedImage->dim[0] == 2)
    {
      for (unsigned int y = ignoreBoundaySize; y < expectedImage->ny - ignoreBoundaySize; ++y)
      {
        for (unsigned int x = ignoreBoundaySize; x < expectedImage->nx - ignoreBoundaySize; ++x)
        {
          int idx = y * expectedImage->nx + x;
          VoxelType absDif = std::abs( imageDataPter[idx] - expectedImageDataPter[idx] );

          if (absDif > allowedVoxelDeviation)
          {
            // Make sure that if both images are not a number, then do something else
            if (imageDataPter[idx] != imageDataPter[idx] || expectedImageDataPter[idx] != expectedImageDataPter[idx]) continue;

            maxDeviation = maxDeviation < absDif ? absDif : maxDeviation;
            ++numberOfDeviations;
          }
        }
      }
    }
  }
  else
  {
    for (unsigned int i = 0; i < expectedImage->nvox; ++i)
    {
      VoxelType absDif = std::abs( imageDataPter[i] - expectedImageDataPter[i] );
      if (absDif > allowedVoxelDeviation)
      {
        // Make sure that if both images are not a number, then do something else
        if (imageDataPter[i] != imageDataPter[i] || expectedImageDataPter[i] != expectedImageDataPter[i]) continue;

        maxDeviation = maxDeviation < absDif ? absDif : maxDeviation;
        ++numberOfDeviations;
      }
    }
  }
  if ( numberOfDeviations > maxAllowedNumberOfDeviations)
  {
    imagesEqual = false;
    char error_msg[100];
    sprintf( error_msg, "Image contents not equal, %i differences larger than %f\n", numberOfDeviations, allowedVoxelDeviation );
    sprintf( error_msg, "Maximum difference between images is %f", maxDeviation );
    supremo_print_error( error_msg );
  }

  // Return the correct return value
  if ( imagesEqual )
  {
    return EXIT_SUCCESS;
  }
  else return EXIT_FAILURE;
}




/** Test if two float values are almost equal. 
 *  First checks if the absolute difference is within the given maxDiff range. If this
 *  is not the case, the relative difference will be checked. 
 * 
 *  \param A first value, compared to second
 *  \param B second value, compared to first
 *  \param maxDiff The maximum absolute difference allowed between values A and B. 
 *  \param maxRelDiff The maximum relative difference between values A and B
 */
bool AlmostEqualRelativeAndAbs( float A, float B,
  float maxDiff, float maxRelDiff = FLT_EPSILON )
{
  // Check if the numbers are really close -- needed
  // when comparing numbers near zero.
  float diff = fabs( A - B );
  if (diff <= maxDiff)
    return true;

  A = fabs( A );
  B = fabs( B );
  float largest = (B > A) ? B : A;

  if (diff <= largest * maxRelDiff)
    return true;
  return false;
}
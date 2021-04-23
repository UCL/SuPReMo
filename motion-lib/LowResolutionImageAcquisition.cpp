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





#include "LowResolutionImageAcquisition.h"
#include <cmath>


//--------------------------------------------------------------
// LowResolutionImageAcquisition::LowResolutionImageAcquisition
//--------------------------------------------------------------
LowResolutionImageAcquisition::LowResolutionImageAcquisition() : 
  lowResolutionThreshold( 1.1f ), 
  roundErrorThreshold( 0.001f ),
  paddingValue( std::numeric_limits<Transformation::PrecisionType>::quiet_NaN() )
{}




//---------------------------------------------------------------
// LowResolutionImageAcquisition::~LowResolutionImageAcquisition
//---------------------------------------------------------------
LowResolutionImageAcquisition::~LowResolutionImageAcquisition()
{}




//---------------------------------------------------------
// LowResolutionImageAcquisition::SimulateImageAcquisition
//---------------------------------------------------------
nifti_image* LowResolutionImageAcquisition::SimulateImageAcquisition( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace )
{

  // Allocate the simualted dynamic image (or low resolution image), which must have the same size as the image in acquisition space
  nifti_image* lowResImage = nifti_copy_nim_info( imgInAcquisitionSpace );
  lowResImage->data = (void *)calloc( lowResImage->nvox, lowResImage->nbyper );

  // Define some pseudonyms for simplicity
  nifti_image* highResImage = imgInFullImgSpace;

  // Check the input data for data-type and dimensinality > 3
  if (lowResImage->datatype != highResImage->datatype)
  {
    supremo_print_error( "High- and low-resolution images have to be of same data-type" );
    supremo_exit( -1, __FILE__, __LINE__ );
  }  
  if ( lowResImage->nt*lowResImage->nu*lowResImage->nv*lowResImage->nw !=
       highResImage->nt*highResImage->nu*highResImage->nv*highResImage->nw)
  {
    supremo_print_error( "High- and low-resolution images have have the same dimensions > 3." );
    supremo_exit( -1, __FILE__, __LINE__ );
  }
   
  // The high-resolution image is resampled onto lower-resolution image
  // for dimensions where low-res voxel size / high-res voxel size > lowResThresh
  // a Gaussian kernel is used to represent the resampling PSF
  // linear interpolation is used for other dimensions
  bool useGaussian[3];
  float ratio, variance, sigma, gaussConst1[3], gaussConst2[3], kernelWidth[3];
  useGaussian[0] = useGaussian[1] = useGaussian[2] = false;
  kernelWidth[0] = kernelWidth[1] = kernelWidth[2] = 1.0;
  variance = sigma = gaussConst1[0] = gaussConst1[1] = gaussConst1[2] = gaussConst2[0] = gaussConst2[1] = gaussConst2[2] = 0.0;
  
  ratio = lowResImage->dx / highResImage->dx;
  
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[0] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[0] = sigma * 2.506628274631;
    gaussConst2[0] = 2 * variance;
    kernelWidth[0] = 3 * sigma;
    //note - variables are in terms of high res voxels, e.g. to get sigma in mm -> sigma * highResImage->dx
  }
  
  ratio = lowResImage->dy / highResImage->dy;
  
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[1] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[1] = sigma * 2.506628274631;
    gaussConst2[1] = 2 * variance;
    kernelWidth[1] = 3 * sigma;
  }
  
  ratio = lowResImage->dz / highResImage->dz;
  
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[2] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[2] = sigma * 2.506628274631;
    gaussConst2[2] = 2 * variance;
    kernelWidth[2] = 3 * sigma;
  }


  // some useful variables and pointers
  float* lowResPtr = static_cast<float*>(lowResImage->data);
  float* highResPtr = static_cast<float*>( highResImage->data);
  int numLowResVox = lowResImage->nx * lowResImage->ny * lowResImage->nz;
  int numHighResVox = highResImage->nx * highResImage->ny * highResImage->nz;
  int numTUVWDims = lowResImage->nt * lowResImage->nu * lowResImage->nv * lowResImage->nw;
  PrecisionType padVal = this->paddingValue;

  // get spacing and origin of highResImage in lowResImage
  mat44 lowResToHighResMatrix, lowResXYZMatrix, highResIJKMatrix;
  
  if (highResImage->sform_code > 0)
    highResIJKMatrix = highResImage->sto_ijk;
  else
    highResIJKMatrix = highResImage->qto_ijk;
  
  if (lowResImage->sform_code > 0)
    lowResXYZMatrix = lowResImage->sto_xyz;
  else
    lowResXYZMatrix = lowResImage->qto_xyz;

  lowResToHighResMatrix = reg_mat44_mul( &(highResIJKMatrix), &(lowResXYZMatrix) );
  float lowResSpacingInHighRes[3], lowResOriginInHighRes[3];
  lowResSpacingInHighRes[0] = lowResToHighResMatrix.m[0][0];
  lowResSpacingInHighRes[1] = lowResToHighResMatrix.m[1][1];
  lowResSpacingInHighRes[2] = lowResToHighResMatrix.m[2][2];
  lowResOriginInHighRes[0] = lowResToHighResMatrix.m[0][3];
  lowResOriginInHighRes[1] = lowResToHighResMatrix.m[1][3];
  lowResOriginInHighRes[2] = lowResToHighResMatrix.m[2][3];

  // loop over dims > 3
  for (int tuvwInd = 0; tuvwInd < numTUVWDims; tuvwInd++)
  {
    //pointers to image data for this tuvw
    float* lowResTUVWPtr = &lowResPtr[tuvwInd*numLowResVox];
    float* highResTUVWPtr = &highResPtr[tuvwInd*numHighResVox];

    // convolve in z direction
    // need variable to hold temporary results - initialise with 0
    int numTmp1Vox = lowResImage->nz*highResImage->ny*highResImage->nx;
    float *tmpVals1 = new float[numTmp1Vox];
    int indTmp1;
    for (indTmp1 = 0; indTmp1 < numTmp1Vox; indTmp1++)
    {
      tmpVals1[indTmp1] = 0.0;
    }

    // declare variables used in omp loop
    int zLowRes, firstHighRes, lastHighRes, indK, indHighRes, zHighRes, yHighRes, xHighRes;
    float lowResCoordInHighRes, *kernel, kDist, kSum;

    // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, highResImage, useGaussian, gaussConst1, gaussConst2, kernelWidth, padVal, \
	highResTUVWPtr, lowResSpacingInHighRes, lowResOriginInHighRes, tmpVals1) \
	private(zLowRes, firstHighRes, lastHighRes, indK, indHighRes, indTmp1, zHighRes, yHighRes, xHighRes, \
	lowResCoordInHighRes, kernel, kDist, kSum)
#endif

    // loop over z in low res image
    for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
    {
      // find coord of zLowRes in high res image
      lowResCoordInHighRes = (float)zLowRes * lowResSpacingInHighRes[2] + lowResOriginInHighRes[2];
      
      // check for rounding error
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      // calculate kernel values
      //
      // first check high res image is 3D
      if (highResImage->nz == 1)
      {
        firstHighRes = lastHighRes = 0;
        kernel = new float[1];
        kernel[0] = 1.0;
      }
      else
      {
        firstHighRes = int( std::ceil( lowResCoordInHighRes - kernelWidth[2] ) );
        lastHighRes = int( std::floor( lowResCoordInHighRes + kernelWidth[2] ) );
        kernel = new float[lastHighRes - firstHighRes + 1];
        
        // Gaussian or linear?
        if (useGaussian[2])
        {
          // Gaussian
          kSum = 0.0;
          for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
          {
            kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
            kernel[indK] = exp( -(kDist*kDist) / gaussConst2[2] ) / gaussConst1[2];
            kSum += kernel[indK];
          }
          for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
          {
            kernel[indK] /= kSum;
          }
        }
        else
        {
          // linear
          // check if lowResCoordInHighRes is exact voxel centre
          if (lastHighRes - firstHighRes > 1)
          {
            firstHighRes = lastHighRes = firstHighRes + 1;
            delete[] kernel;
            kernel = new float[1];
            kernel[0] = 1.0;
          }
          else
          {
            kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
            kernel[0] = 1.0 - kDist;
            kernel[1] = kDist;
          }
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // index into tmp1
        indTmp1 = zLowRes * highResImage->ny * highResImage->nx;

        // check high res vox in FOV
        zHighRes = firstHighRes + indK;
        if (zHighRes >= 0 && zHighRes < highResImage->nz)
        {
          // index into high res image
          indHighRes = zHighRes * highResImage->ny*highResImage->nx;

          // loop over high res vox and store intermediate results
          for (xHighRes = 0; xHighRes < highResImage->nx; xHighRes++)
          {
            for (yHighRes = 0; yHighRes < highResImage->ny; yHighRes++)
            {
              // if NaN value replace with paddingValue
              if (!isnan( static_cast<float>(highResTUVWPtr[indHighRes]) ))
              {
                tmpVals1[indTmp1] += kernel[indK] * static_cast<float>(highResTUVWPtr[indHighRes]);
              }
              else
              {
                tmpVals1[indTmp1] += kernel[indK] * static_cast<float>(padVal);
              }
              // increment indices
              indTmp1++;
              indHighRes++;
            }
          }// loop over high res vox and store intermediate results
        }// check high res vox in FOV
        else
        {
          // if outside FOV use paddingValue
          // loop over high res vox and store intermediate results
          for (xHighRes = 0; xHighRes < highResImage->nx; xHighRes++)
          {
            for (yHighRes = 0; yHighRes < highResImage->ny; yHighRes++)
            {
              tmpVals1[indTmp1] += kernel[indK] * static_cast<float>(padVal);
              // increment indices
              indTmp1++;
            }
          }// loop over high res vox and store intermediate results
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over z in low res image (OMP loop)


    // convolve in y direction
    // need variable to hold temporary results - initialise with 0
    int numTmp2Vox = lowResImage->nz*lowResImage->ny*highResImage->nx;
    float *tmpVals2 = new float[numTmp2Vox];
    int indTmp2;
    for (indTmp2 = 0; indTmp2 < numTmp2Vox; indTmp2++)
    {
      tmpVals2[indTmp2] = 0.0;
    }

    // declare new variables used in omp loop
    int yLowRes;

    // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, highResImage, useGaussian, gaussConst1, gaussConst2, kernelWidth, padVal, \
	lowResSpacingInHighRes, lowResOriginInHighRes, tmpVals1, tmpVals2) \
	private(zLowRes, yLowRes, firstHighRes, lastHighRes, indK, indTmp1, indTmp2, yHighRes, xHighRes, \
	lowResCoordInHighRes, kernel, kDist, kSum)
#endif

    // loop over y in low res image
    for (yLowRes = 0; yLowRes < lowResImage->ny; yLowRes++)
    {
      // find coord of yLowRes in high res image
      lowResCoordInHighRes = (float)yLowRes*lowResSpacingInHighRes[1] + lowResOriginInHighRes[1];
      // check for rounding error
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      // calculate kernel values
      firstHighRes = int( std::ceil( lowResCoordInHighRes - kernelWidth[1] ) );
      lastHighRes = int( std::floor( lowResCoordInHighRes + kernelWidth[1] ) );
      kernel = new float[lastHighRes - firstHighRes + 1];
      // Gaussian or linear?
      if (useGaussian[1])
      {
        // Gaussian
        kSum = 0.0;
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
          kernel[indK] = exp( -(kDist*kDist) / gaussConst2[1] ) / gaussConst1[1];
          kSum += kernel[indK];
        }
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kernel[indK] /= kSum;
        }
      }
      else
      {
        // linear
        // check if lowResCoordInHighRes is exact voxel centre
        if (lastHighRes - firstHighRes > 1)
        {
          firstHighRes = lastHighRes = firstHighRes + 1;
          delete[] kernel;
          kernel = new float[1];
          kernel[0] = 1.0;
        }
        else
        {
          kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
          kernel[0] = 1.0 - kDist;
          kernel[1] = kDist;
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // check high res vox in FOV
        yHighRes = firstHighRes + indK;
        if (yHighRes >= 0 && yHighRes < highResImage->ny)
        {
          // loop over values in tmpVals1 and store intermediate results in tmpVals2
          for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
          {
            // indices into tmpVals1 and tmpVals2
            indTmp1 = (zLowRes * highResImage->ny + yHighRes) * highResImage->nx;
            indTmp2 = (zLowRes * lowResImage->ny + yLowRes) * highResImage->nx;

            for (xHighRes = 0; xHighRes < highResImage->nx; xHighRes++)
            {
              // no need to test for NaN values anymore as they will only be in tmpVals1
              // if paddingValue = NaN
              tmpVals2[indTmp2] += kernel[indK] * tmpVals1[indTmp1];
              indTmp1++;
              indTmp2++;
            }
          }
        }// check high res vox in FOV
        else
        {
          // if outside FOV use paddingValue
          // loop over values in tmpVals1 and store intermediate results in tmpVals2
          for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
          {
            //indices into tmpVals1 and tmpVals2
            indTmp2 = (zLowRes * lowResImage->ny + yLowRes) * highResImage->nx;

            for (xHighRes = 0; xHighRes < highResImage->nx; xHighRes++)
            {
              // no need to test for NaN values anymore as they will only be in tmpVals1
              // if paddingValue = NaN
              tmpVals2[indTmp2] += kernel[indK] * static_cast<float>(padVal);
              indTmp2++;
            }
          }
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over y in low res image (OMP loop)

    // finished with tmpVals1 so delete
    delete[] tmpVals1;
    tmpVals1 = NULL;


    // convolve in x direction
    // need float variable to hold results, will be converted to ImageTYPE at end
    // initialise with 0
    float *lowResData = new float[numLowResVox];
    int indLowRes;
    for (indLowRes = 0; indLowRes < numLowResVox; indLowRes++)
    {
      lowResData[indLowRes] = 0.0;
    }

    // declare new variables used in omp loop
    int xLowRes;

    // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, highResImage, useGaussian, gaussConst1, gaussConst2, kernelWidth, padVal, \
	lowResSpacingInHighRes, lowResOriginInHighRes, tmpVals2, lowResData) \
	private(zLowRes, yLowRes, xLowRes, firstHighRes, lastHighRes, indK, indTmp2, indLowRes, xHighRes, \
	lowResCoordInHighRes, kernel, kDist, kSum)
#endif

    // loop over x in low res image
    for (xLowRes = 0; xLowRes < lowResImage->nx; xLowRes++)
    {
      // find coord of xLowRes in high res image
      lowResCoordInHighRes = (float)xLowRes*lowResSpacingInHighRes[0] + lowResOriginInHighRes[0];
      // check for rounding error
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      // calculate kernel values
      firstHighRes = int( std::ceil( lowResCoordInHighRes - kernelWidth[0] ) );
      lastHighRes = int( std::floor( lowResCoordInHighRes + kernelWidth[0] ) );
      kernel = new float[lastHighRes - firstHighRes + 1];
      // Gaussian or linear?
      if (useGaussian[0])
      {
        // Gaussian
        kSum = 0.0;
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
          kernel[indK] = exp( -(kDist*kDist) / gaussConst2[0] ) / gaussConst1[0];
          kSum += kernel[indK];
        }
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kernel[indK] /= kSum;
        }
      }
      else
      {
        // linear
        // check if lowResCoordInHighRes is exact voxel centre
        if (lastHighRes - firstHighRes > 1)
        {
          firstHighRes = lastHighRes = firstHighRes + 1;
          delete[] kernel;
          kernel = new float[1];
          kernel[0] = 1.0;
        }
        else
        {
          kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
          kernel[0] = 1.0 - kDist;
          kernel[1] = kDist;
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // check high res vox in FOV
        xHighRes = firstHighRes + indK;
        if (xHighRes >= 0 && xHighRes < highResImage->nx)
        {
          // loop over values in tmpVals2 and store results in lowResData
          for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
          {
            for (yLowRes = 0; yLowRes < lowResImage->ny; yLowRes++)
            {
              indTmp2 = (zLowRes * lowResImage->ny + yLowRes) * highResImage->nx + xHighRes;
              indLowRes = (zLowRes * lowResImage->ny + yLowRes) * lowResImage->nx + xLowRes;
              lowResData[indLowRes] += kernel[indK] * tmpVals2[indTmp2];
            }
          }
        }// check high res vox in FOV
        else
        {
          // if outside FOV use paddingValue
          // loop over values in tmpVals2 and store results in lowResData
          for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
          {
            for (yLowRes = 0; yLowRes < lowResImage->ny; yLowRes++)
            {
              indLowRes = (zLowRes * lowResImage->ny + yLowRes) * lowResImage->nx + xLowRes;
              lowResData[indLowRes] += kernel[indK] * static_cast<float>(padVal);
            }
          }
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over x in low res image (OMP loop)

    // finished with tmpVals2 so delete
    delete[] tmpVals2;
    tmpVals2 = NULL;

    // copy lowResData into low res image converting to ImageTYPE
    for (indLowRes = 0; indLowRes < numLowResVox; indLowRes++)
    {
      lowResTUVWPtr[indLowRes] = (float)lowResData[indLowRes];
    }

    //finished with lowResData so delete
    delete[] lowResData;
    lowResData = NULL;

  }// loop over dims > 3
  
  return lowResImage;
}




//-------------------------------------------------
// LowResolutionImageAcquisition::CalculateAdjoint
//-------------------------------------------------
void LowResolutionImageAcquisition::CalculateAdjoint( nifti_image* imgInFullImgSpace, nifti_image* imgInAcquisitionSpace )
{
  // The current image in acquisition space is remembered for the allocation of the image after adjoint
  this->curImageInAcquisitionSpace = imgInAcquisitionSpace;
  this->curImageInFullImgSpace = imgInFullImgSpace;

  // Allocate the image after applying the adjoint and the corresponding weights image
  this->AllocateImageAfterAdjoint();
  this->AllocateWeightsImageAfterAdjoint();


  // Intermediate translation map        
  // lowResImage  == currentReference    |== this->curImageInAcquisitionSpace
  // highResImage == dynamicAfterAdjoint |== this->imageAfterAdjoint
  // lowResMask   == currentMask         |== N/A
  // kernelSumMap == weightsInDynImg     |== this->weightsImageAfterAdjoint

  // For simplicity rename the pointer to the images. 
  nifti_image* lowResImage = this->curImageInAcquisitionSpace;
  nifti_image* highResImage = this->imageAfterAdjoint;

  // the low-res image is 'spread out' onto the high res image
  // for dimensions where low-res voxel size / high-res voxel size > lowResThresh
  // a gaussian kernel is used to represent the resampling PSF
  // linear interpolation is used for other dimensions
  bool useGaussian[3];
  float ratio, variance, sigma, gaussConst1[3], gaussConst2[3], kernelWidth[3];
  useGaussian[0] = useGaussian[1] = useGaussian[2] = false;
  kernelWidth[0] = kernelWidth[1] = kernelWidth[2] = 1.0;
  variance = sigma = gaussConst1[0] = gaussConst1[1] = gaussConst1[2] = gaussConst2[0] = gaussConst2[1] = gaussConst2[2] = 0.0;

  ratio = lowResImage->dx / highResImage->dx;
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[0] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[0] = sigma * 2.506628274631;
    gaussConst2[0] = 2 * variance;
    kernelWidth[0] = 3 * sigma;
    //note - variables are in terms of high res voxels, e.g. to get sigma in mm -> sigma * highResImage->dx
  }

  ratio = lowResImage->dy / highResImage->dy;
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[1] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[1] = sigma * 2.506628274631;
    gaussConst2[1] = 2 * variance;
    kernelWidth[1] = 3 * sigma;
  }

  ratio = lowResImage->dz / highResImage->dz;
  if (ratio > this->lowResolutionThreshold)
  {
    useGaussian[2] = true;
    variance = (std::pow( ratio, 2 ) - 1) * 0.18033688; // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    sigma = sqrt( variance );
    gaussConst1[2] = sigma * 2.506628274631;
    gaussConst2[2] = 2 * variance;
    kernelWidth[2] = 3 * sigma;
  }

  // some useful variables and pointers
  float *lowResPtr = static_cast<float*>(lowResImage->data);
  float *highResPtr = static_cast<float*>(highResImage->data);
  int numLowResVox = lowResImage->nx * lowResImage->ny * lowResImage->nz;
  int numHighResVox = highResImage->nx * highResImage->ny * highResImage->nz;
  int numTUVWDims = lowResImage->nt * lowResImage->nu * lowResImage->nv * lowResImage->nw;

  // get spacing and origin of highResImage in lowResImage
  mat44 lowResToHighResMatrix, lowResXYZMatrix, highResIJKMatrix;
  if (highResImage->sform_code > 0)
    highResIJKMatrix = highResImage->sto_ijk;
  else
    highResIJKMatrix = highResImage->qto_ijk;

  if (lowResImage->sform_code > 0)
    lowResXYZMatrix = lowResImage->sto_xyz;
  else
    lowResXYZMatrix = lowResImage->qto_xyz;

  // Generate a matrix that calculates corresponding indices from the low res image to the high res image
  // Essentially combine the two multiplications 
  // 1) lowResXYZMatrix * (i,j,k)_low = (x,y,z)  (physical coordinate centered on low resolution voxel)
  // 2) highResIJKMatrix * (x,y,z) = (i,j,k)_high
  // into one
  // (highResIJKMatrix * lowResXYZMatrix) * (i,j,k)_low = (i,j,k)_high
  lowResToHighResMatrix = reg_mat44_mul( &(highResIJKMatrix), &(lowResXYZMatrix) );

  // Note: The calculation of the spacing between low-resolution and high resolution images by only considering the main diagonal 
  //       requires no permutation of directions between those images.
  float lowResSpacingInHighRes[3], lowResOriginInHighRes[3];
  lowResSpacingInHighRes[0] = lowResToHighResMatrix.m[0][0];
  lowResSpacingInHighRes[1] = lowResToHighResMatrix.m[1][1];
  lowResSpacingInHighRes[2] = lowResToHighResMatrix.m[2][2];
  lowResOriginInHighRes[0] = lowResToHighResMatrix.m[0][3];
  lowResOriginInHighRes[1] = lowResToHighResMatrix.m[1][3];
  lowResOriginInHighRes[2] = lowResToHighResMatrix.m[2][3];

  float* kernelSumMap = static_cast<float*>(this->weightsImageAfterAdjoint->data);

  // loop over dims > 3
  for (int tuvwInd = 0; tuvwInd < numTUVWDims; tuvwInd++)
  {
    // pointers to image data for this tuvw
    float *lowResTUVWPtr = &lowResPtr[tuvwInd*numLowResVox];
    float *highResTUVWPtr = &highResPtr[tuvwInd*numHighResVox];

    // convolve in z direction
    // need variables to hold temporary results - initialise with 0
    int numTmp1Vox = highResImage->nz*lowResImage->ny*lowResImage->nx;
    float *tmpVals1 = new float[numTmp1Vox];
    float *tmpKSM1 = new float[numTmp1Vox];
    int indTmp1;
    for (indTmp1 = 0; indTmp1 < numTmp1Vox; indTmp1++)
    {
      tmpVals1[indTmp1] = 0.0;
      tmpKSM1[indTmp1] = 0.0;
    }

    // declare variables used in omp loop
    int zLowRes, yLowRes, xLowRes, firstHighRes, lastHighRes, indK, indLowRes, zHighRes;
    float lowResCoordInHighRes, *kernel, kDist, kSum;

    // loop over z in low res image
    for (zLowRes = 0; zLowRes < lowResImage->nz; zLowRes++)
    {
      // find coord of zLowRes in high res image
      lowResCoordInHighRes = (float)zLowRes * lowResSpacingInHighRes[2] + lowResOriginInHighRes[2];

      // check if we can safely round without losing too much precision
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      //calculate kernel values
      //
      //first check high res image is 3D
      if (highResImage->nz == 1)
      {
        firstHighRes = lastHighRes = 0;
        kernel = new float[1];
        kernel[0] = 1.0;
      }
      else
      {
        firstHighRes = (int)(std::ceil( lowResCoordInHighRes - kernelWidth[2] ));
        lastHighRes = (int)(std::floor( lowResCoordInHighRes + kernelWidth[2] ));
        kernel = new float[lastHighRes - firstHighRes + 1];
        // Gaussian or linear?
        if (useGaussian[2])
        {
          // Gaussian
          kSum = 0.0;
          for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
          {
            kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
            kernel[indK] = exp( -(kDist*kDist) / gaussConst2[2] ) / gaussConst1[2];
            kSum += kernel[indK];
          }
          for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
          {
            kernel[indK] /= kSum;
          }
        }
        else
        {
          // linear

          // check if lowResCoordInHighRes is exact voxel centre
          if (lastHighRes - firstHighRes > 1)
          {
            firstHighRes = lastHighRes = firstHighRes + 1;
            delete[] kernel;
            kernel = new float[1];
            kernel[0] = 1.0;
          }
          else
          {
            kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
            kernel[0] = 1.0 - kDist;
            kernel[1] = kDist;
          }
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // check high res vox in FOV (if not no need to do anything)
        zHighRes = firstHighRes + indK;
        if (zHighRes >= 0 && zHighRes < highResImage->nz)
        {
          // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, tmpVals1, tmpKSM1, zLowRes, zHighRes, indK, kernel, lowResTUVWPtr) \
	private(yLowRes, xLowRes, indLowRes, indTmp1)
#endif
          // loop over low res vox and store intermediate results
          for (yLowRes = 0; yLowRes < lowResImage->ny; yLowRes++)
          {
            // indices into lowResImage and tmp1
            indLowRes = (zLowRes * lowResImage->ny + yLowRes) * lowResImage->nx;
            indTmp1 = (zHighRes * lowResImage->ny + yLowRes) * lowResImage->nx;

            for (xLowRes = 0; xLowRes < lowResImage->nx; xLowRes++)
            {
              // not NaN
              if (!isnan( static_cast<float>(lowResTUVWPtr[indLowRes]) ))
              {
                tmpVals1[indTmp1] += kernel[indK] * static_cast<float>(lowResTUVWPtr[indLowRes]);
                tmpKSM1[indTmp1] += kernel[indK];
              }
              // increment indices
              indTmp1++;
              indLowRes++;
            }
          }// loop over low res vox and store intermediate results (OMP loop)
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over z in low res image

    // convolve in y direction

    // need variables to hold temporary results - initialise with 0
    int numTmp2Vox = highResImage->nz*highResImage->ny*lowResImage->nx;
    float *tmpVals2 = new float[numTmp2Vox];
    float *tmpKSM2 = new float[numTmp2Vox];
    int indTmp2;
    
    for (indTmp2 = 0; indTmp2 < numTmp2Vox; indTmp2++)
    {
      tmpVals2[indTmp2] = 0.0;
      tmpKSM2[indTmp2] = 0.0;
    }

    // declare new variables used in omp loop
    int yHighRes;

    // loop over y in low res image
    for (yLowRes = 0; yLowRes < lowResImage->ny; yLowRes++)
    {
      // find coord of yLowRes in high res image
      lowResCoordInHighRes = (float)yLowRes * lowResSpacingInHighRes[1] + lowResOriginInHighRes[1];
      // check if we can round without losing precision
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      // calculate kernel values
      firstHighRes = (int)( std::ceil( lowResCoordInHighRes - kernelWidth[1] ) );
      lastHighRes = (int)( std::floor( lowResCoordInHighRes + kernelWidth[1] ) );
      kernel = new float[lastHighRes - firstHighRes + 1];
      // Gaussian or linear?
      if (useGaussian[1])
      {
        // Gaussian
        kSum = 0.0;
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
          kernel[indK] = exp( -(kDist*kDist) / gaussConst2[1] ) / gaussConst1[1];
          kSum += kernel[indK];
        }
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kernel[indK] /= kSum;
        }
      }
      else
      {
        // linear
        // check if lowResCoordInHighRes is exact voxel centre
        if (lastHighRes - firstHighRes > 1)
        {
          firstHighRes = lastHighRes = firstHighRes + 1;
          delete[] kernel;
          kernel = new float[1];
          kernel[0] = 1.0;
        }
        else
        {
          kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
          kernel[0] = 1.0 - kDist;
          kernel[1] = kDist;
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // check high res vox in FOV (if not no need to do anything)
        yHighRes = firstHighRes + indK;
        if (yHighRes >= 0 && yHighRes < highResImage->ny)
        {
          // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, highResImage, tmpVals1, tmpKSM1, tmpVals2, tmpKSM2, yLowRes, yHighRes, indK, kernel) \
	private(zHighRes, xLowRes, indTmp1, indTmp2)
#endif
          // loop over values in tmpVals1 and store intermediate results in tmpVals2
          for (zHighRes = 0; zHighRes < highResImage->nz; zHighRes++)
          {
            // indices into tmpVals1 and tmpVals2
            indTmp1 = (zHighRes * lowResImage->ny + yLowRes) * lowResImage->nx;
            indTmp2 = (zHighRes * highResImage->ny + yHighRes) * lowResImage->nx;

            for (xLowRes = 0; xLowRes < lowResImage->nx; xLowRes++)
            {
              tmpVals2[indTmp2] += kernel[indK] * tmpVals1[indTmp1];
              tmpKSM2[indTmp2] += kernel[indK] * tmpKSM1[indTmp1];
              indTmp1++;
              indTmp2++;
            }
          }// loop over values in tmpVals1 and store intermediate results in tmpVals2 (OMP loop)
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over y in low res image

    // finished with tmpVals1 and tmpKSM1 so delete
    delete[] tmpVals1;
    tmpVals1 = NULL;
    delete[] tmpKSM1;
    tmpKSM1 = NULL;


    // convolve in x direction

    // need float variable to hold results, will be converted to ImageTYPE at end
    // initialise with 0
    // also initialise kernelSumMap with 0
    float *highResData = new float[numHighResVox];
    int indHighRes;
    for (indHighRes = 0; indHighRes < numHighResVox; indHighRes++)
    {
      highResData[indHighRes] = 0.0;
      kernelSumMap[indHighRes] = 0.0;
    }

    // declare new variables used in omp loop
    int xHighRes;

    // loop over x in low res image
    for (xLowRes = 0; xLowRes < lowResImage->nx; xLowRes++)
    {
      // find coord of xLowRes in high res image
      lowResCoordInHighRes = (float)xLowRes*lowResSpacingInHighRes[0] + lowResOriginInHighRes[0];
      // check if we can round without losing accuracy
      if (fabs( lowResCoordInHighRes - std::round( lowResCoordInHighRes ) ) < this->roundErrorThreshold)
      {
        lowResCoordInHighRes = std::round( lowResCoordInHighRes );
      }

      // calculate kernel values
      firstHighRes = int( std::ceil( lowResCoordInHighRes - kernelWidth[0] ) );
      lastHighRes = int( std::floor( lowResCoordInHighRes + kernelWidth[0] ) );
      kernel = new float[lastHighRes - firstHighRes + 1];
      // Gaussian or linear?
      if (useGaussian[0])
      {
        // Gaussian
        kSum = 0.0;
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kDist = static_cast<float>(firstHighRes + indK) - lowResCoordInHighRes;
          kernel[indK] = exp( -(kDist*kDist) / gaussConst2[0] ) / gaussConst1[0];
          kSum += kernel[indK];
        }
        for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
        {
          kernel[indK] /= kSum;
        }
      }
      else
      {
        // linear
        // check if lowResCoordInHighRes is exact voxel centre
        if (lastHighRes - firstHighRes > 1)
        {
          firstHighRes = lastHighRes = firstHighRes + 1;
          delete[] kernel;
          kernel = new float[1];
          kernel[0] = 1.0;
        }
        else
        {
          kDist = lowResCoordInHighRes - static_cast<float>(firstHighRes);
          kernel[0] = 1.0 - kDist;
          kernel[1] = kDist;
        }
      }// calculate kernel values

      // loop over kernel
      for (indK = 0; indK <= lastHighRes - firstHighRes; indK++)
      {
        // check high res vox in FOV (if not no need to do anything)
        xHighRes = firstHighRes + indK;
        if (xHighRes >= 0 && xHighRes < highResImage->nx)
        {
          // use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(lowResImage, highResImage, tmpVals2, tmpKSM2, highResData, kernelSumMap, xLowRes, xHighRes, indK, kernel) \
	private(zHighRes, yHighRes, indTmp2, indHighRes)
#endif					
          // loop over values in tmpVals2 and store results in highResData
          for (zHighRes = 0; zHighRes < highResImage->nz; zHighRes++)
          {
            for (yHighRes = 0; yHighRes < highResImage->ny; yHighRes++)
            {
              indTmp2 = (zHighRes * highResImage->ny + yHighRes) * lowResImage->nx + xLowRes;
              indHighRes = (zHighRes * highResImage->ny + yHighRes) * highResImage->nx + xHighRes;
              highResData[indHighRes] += kernel[indK] * tmpVals2[indTmp2];
              kernelSumMap[indHighRes] += kernel[indK] * tmpKSM2[indTmp2];
            }
          }// loop over values in tmpVals2 and store results in highResData (OMP loop)
        }// check high res vox in FOV
      }// loop over kernel
      delete[] kernel;
    }// loop over x in low res image

    // finished with tmpVals2 and tmpKSM2 so delete
    delete[] tmpVals2;
    tmpVals2 = NULL;
    delete[] tmpKSM2;
    tmpKSM2 = NULL;

    // copy highResData into high res image
    for (indHighRes = 0; indHighRes < numHighResVox; indHighRes++)
    {
      highResTUVWPtr[indHighRes] = (float) highResData[indHighRes];
    }

    // finished with highResData so delete
    delete[] highResData;
    highResData = NULL;
  }// loop over dims > 3
}



//---------------------------------------------------------------------
// LowResolutionImageAcquisition::AllocateMinimumSizeImgInFullImgSpace
//---------------------------------------------------------------------
nifti_image * LowResolutionImageAcquisition::AllocateMinimumSizeImgInFullImgSpace( nifti_image * imgInFullImgSpace, nifti_image * imgInAcquisitionSpace )
{
  // The size of the allocated image needs to be determined. Starting with size is the image in acquisition space (the lower resolution),
  // the size will be refined below using the geometry information provided by imgInFullImgSpace (the higher resolution).
  // Memory allocation will happen only after size was determined.
  
  // Allocate the image
  nifti_image* minSizedImgInFullImgSpace;

  // Fill in the header information 
  minSizedImgInFullImgSpace = nifti_copy_nim_info( imgInAcquisitionSpace );
  
  // The offset in voxels in every image dimension. Required if the image
  float originOffsetVox[3];
  originOffsetVox[0] = originOffsetVox[1] = originOffsetVox[2] = 0.0;

  // Determin image size in x-dimension
  if (minSizedImgInFullImgSpace->dx > imgInFullImgSpace->dx * this->lowResolutionThreshold)
  {
    // Update voxel size 
    minSizedImgInFullImgSpace->dx = minSizedImgInFullImgSpace->pixdim[1] = imgInFullImgSpace->dx;

    // Also voxel size in sform
    minSizedImgInFullImgSpace->sto_xyz.m[0][0] = imgInFullImgSpace->sto_xyz.m[0][0];
    minSizedImgInFullImgSpace->sto_xyz.m[0][1] = imgInFullImgSpace->sto_xyz.m[0][1];
    minSizedImgInFullImgSpace->sto_xyz.m[0][2] = imgInFullImgSpace->sto_xyz.m[0][2];
    
    // Update number of slices to cover sufficient extent
    float sigma_mm = sqrt( (std::pow( imgInAcquisitionSpace->dx, 2 ) - std::pow( imgInFullImgSpace->dx, 2 )) * 0.18033688 ); // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    float extent_low_mm = 3.f * sigma_mm;
    float extent_high_mm = 3.f * sigma_mm + imgInAcquisitionSpace->dx * (float)(imgInAcquisitionSpace->nx - 1);

    //origin offest is always less/equal-to zero, hence we need to subtract the added offset when calculating the complete extent in voxels
    originOffsetVox[0] = -std::ceil( extent_low_mm / imgInFullImgSpace->dx );
    minSizedImgInFullImgSpace->nx = minSizedImgInFullImgSpace->dim[1] = -originOffsetVox[0] + std::ceil( extent_high_mm / imgInFullImgSpace->dx ) + 1;
  }
  
  // Determin image size in z-dimension
  if (minSizedImgInFullImgSpace->dy > imgInFullImgSpace->dy * this->lowResolutionThreshold)
  {
    // Update voxel size 
    minSizedImgInFullImgSpace->dy = minSizedImgInFullImgSpace->pixdim[2] = imgInFullImgSpace->dy;
    
    // Also voxel size in sform
    minSizedImgInFullImgSpace->sto_xyz.m[1][0] = imgInFullImgSpace->sto_xyz.m[1][0];
    minSizedImgInFullImgSpace->sto_xyz.m[1][1] = imgInFullImgSpace->sto_xyz.m[1][1];
    minSizedImgInFullImgSpace->sto_xyz.m[1][2] = imgInFullImgSpace->sto_xyz.m[1][2];
    
    // Update number of slices to cover sufficient extent
    float sigma_mm = sqrt( (std::pow( imgInAcquisitionSpace->dy, 2 ) - std::pow( imgInFullImgSpace->dy, 2 )) * 0.18033688 ); // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    float extent_low_mm = 3.f * sigma_mm;
    float extent_high_mm = 3.f * sigma_mm + imgInAcquisitionSpace->dy * (float)(imgInAcquisitionSpace->ny - 1);

    //origin offest is always less/equal-to zero, hence we need to subtract the added offset when calculating the complete extent in voxels
    originOffsetVox[1] = -std::ceil( extent_low_mm / imgInFullImgSpace->dy );
    minSizedImgInFullImgSpace->ny = minSizedImgInFullImgSpace->dim[2] = -originOffsetVox[1] + std::ceil( extent_high_mm / imgInFullImgSpace->dy ) + 1;
  }

  // Determin image size in z-dimension
  if (imgInFullImgSpace->nz > 1 && minSizedImgInFullImgSpace->dz > imgInFullImgSpace->dz * this->lowResolutionThreshold)
  {
    // Update voxel size 
    minSizedImgInFullImgSpace->dz = minSizedImgInFullImgSpace->pixdim[3] = imgInFullImgSpace->dz;

    // Also voxel size in sform
    minSizedImgInFullImgSpace->sto_xyz.m[2][0] = imgInFullImgSpace->sto_xyz.m[2][0];
    minSizedImgInFullImgSpace->sto_xyz.m[2][1] = imgInFullImgSpace->sto_xyz.m[2][1];
    minSizedImgInFullImgSpace->sto_xyz.m[2][2] = imgInFullImgSpace->sto_xyz.m[2][2];

    // Update number of slices to cover sufficient extent
    float sigma_mm = sqrt( (std::pow( imgInAcquisitionSpace->dz, 2 ) - std::pow( imgInFullImgSpace->dz, 2 )) * 0.18033688 ); // 0.18033688 = (2*sqrt(2*ln(2)))^-2
    float extent_low_mm = 3.f * sigma_mm;
    float extent_high_mm = 3.f * sigma_mm + imgInAcquisitionSpace->dz * (float)(imgInAcquisitionSpace->nz - 1);

    //origin offest is always less/equal-to zero, hence we need to subtract the added offset when calculating the complete extent in voxels
    originOffsetVox[2] = -std::ceil( extent_low_mm / imgInFullImgSpace->dz );
    minSizedImgInFullImgSpace->nz = minSizedImgInFullImgSpace->dim[3] = -originOffsetVox[2] + std::ceil( extent_high_mm / imgInFullImgSpace->dz ) + 1;
  }

  // Update origin in sform and qform:
  // Convert the new starting index (originOffsetVox) into a physical coordinate
  float originOffsetReal[3];
  reg_mat44_mul( &(minSizedImgInFullImgSpace->sto_xyz), originOffsetVox, originOffsetReal );
  minSizedImgInFullImgSpace->sto_xyz.m[0][3] = minSizedImgInFullImgSpace->qoffset_x = originOffsetReal[0];
  minSizedImgInFullImgSpace->sto_xyz.m[1][3] = minSizedImgInFullImgSpace->qoffset_y = originOffsetReal[1];
  minSizedImgInFullImgSpace->sto_xyz.m[2][3] = minSizedImgInFullImgSpace->qoffset_z = originOffsetReal[2];
  minSizedImgInFullImgSpace->sto_ijk = nifti_mat44_inverse( minSizedImgInFullImgSpace->sto_xyz );
  minSizedImgInFullImgSpace->qto_xyz = nifti_quatern_to_mat44( minSizedImgInFullImgSpace->quatern_b,
    minSizedImgInFullImgSpace->quatern_c,
    minSizedImgInFullImgSpace->quatern_d,
    minSizedImgInFullImgSpace->qoffset_x,
    minSizedImgInFullImgSpace->qoffset_y,
    minSizedImgInFullImgSpace->qoffset_z,
    minSizedImgInFullImgSpace->dx,
    minSizedImgInFullImgSpace->dy,
    minSizedImgInFullImgSpace->dz,
    minSizedImgInFullImgSpace->qfac );
  minSizedImgInFullImgSpace->qto_ijk = nifti_mat44_inverse( minSizedImgInFullImgSpace->qto_xyz );
  
  // allocate data after evaluating the final number of voxels
  minSizedImgInFullImgSpace->nvox = minSizedImgInFullImgSpace->nx
    * minSizedImgInFullImgSpace->ny
    * minSizedImgInFullImgSpace->nz
    * minSizedImgInFullImgSpace->nt
    * minSizedImgInFullImgSpace->nu
    * minSizedImgInFullImgSpace->nv
    * minSizedImgInFullImgSpace->nw;
  minSizedImgInFullImgSpace->data = (void *)calloc( minSizedImgInFullImgSpace->nvox, minSizedImgInFullImgSpace->nbyper );

  return minSizedImgInFullImgSpace;
}




//----------------------------------------------------------
// LowResolutionImageAcquisition::AllocateImageAfterAdjoint
//----------------------------------------------------------
void LowResolutionImageAcquisition::AllocateImageAfterAdjoint()
{
  // Delete the previous image after applying the adjoint of the image acquisition
  this->ClearImageAfterAdjoint();
  this->imageAfterAdjoint = this->AllocateMinimumSizeImgInFullImgSpace(this->curImageInFullImgSpace, this->curImageInAcquisitionSpace);
}


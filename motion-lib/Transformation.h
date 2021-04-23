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


#include "nifti1_io.h"
#include <memory>
#include <vector>


template <class PrecisionType>
void interpolantCubicSpline( PrecisionType ratio, PrecisionType *basis )
{
  if (ratio < 0.0) ratio = 0.0; //reg_rounding error
  PrecisionType FF = ratio * ratio;
  basis[0] = (PrecisionType)((ratio * ((2.0 - ratio)*ratio - 1.0)) / 2.0);
  basis[1] = (PrecisionType)((FF * (3.0*ratio - 5.0) + 2.0) / 2.0);
  basis[2] = (PrecisionType)((ratio * ((4.0 - 3.0*ratio)*ratio + 1.0)) / 2.0);
  basis[3] = (PrecisionType)((ratio - 1.0) * FF / 2.0);
}

template <class PrecisionType>
void interpolantCubicSpline( PrecisionType ratio, PrecisionType *basis, PrecisionType *derivative )
{
  interpolantCubicSpline<PrecisionType>( ratio, basis );
  if (ratio < 0.0) ratio = 0.0; //reg_rounding error
  PrecisionType FF = ratio * ratio;
  derivative[0] = (PrecisionType)((4.0*ratio - 3.0*FF - 1.0) / 2.0);
  derivative[1] = (PrecisionType)((9.0*ratio - 10.0) * ratio / 2.0);
  derivative[2] = (PrecisionType)((8.0*ratio - 9.0*FF + 1) / 2.0);
  derivative[3] = (PrecisionType)((3.0*ratio - 2.0) * ratio / 2.0);
}


/** Transformation class
 *  Abstract base class that must be reimplemented to obtain a functioning transformation class.
 */
class Transformation
{
public:
  typedef float PrecisionType;

  virtual ~Transformation(){};

  /** Initialise the transformation and all internal parameters. 
   *  \param referenceImageIn The reference image defining the space where the transformation will be defined.
   */
  virtual void InitialiseLevel(unsigned int level) = 0;
  
  /** Get a copy of the parameters. The memory will be allocated with malloc and filled by the derived class.
   */
  virtual PrecisionType* GetCopyOfParameters() = 0;

  /** Get the number of parameters required to describe this transformation completely.
   */
  virtual unsigned int GetNumberOfParameters() { return this->numberOfParameters; };

  /** Set the transformation parameters
   */
  virtual void SetParameters( PrecisionType* paramsIn, bool parametersAreDisplacements ) = 0;

  /** Get the deformation vector field according to the current parameters. 
   *  This function is expected to return a pointer to the member variable of the DVF. 
   *  \param targetImageIn Nifti image is used to allocate the DVF. The DVF will remain 
   *                       part of the transformation class and will be deleted if the 
   *                       transformation is deleted.
   */
  virtual nifti_image* GetDeformationVectorField( nifti_image* targetImageIn ) = 0;
  
  /** Compute the gradient of the DVF with respect to the transformation parameters
   *  \f$ \frac{\partial \textbf{DVF}_t}{\partial \mathbf{M_t}} \f$. 
   *  \param denseDVFIn pointer to the nifti image for which the best match of transformation 
   *                    parameters will be found
   *  \return pointer to the internally allocated output data of length numberOfParameters. Needs to be freed outside 
   *          of this class since no internal references will be kept.
   */
  virtual PrecisionType* GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn ) = 0;
  
  /** Get the gradient of the regularisation/constraint with respect to the transformation parameters
   *  \f$ \frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \f$
   */
  virtual PrecisionType* GetConstraintGradientWRTTransformationParameters() = 0;
  
  /** Get the value of the constraint term
   */
  virtual double GetConstraintValue() = 0;
  
  /** Transform/warp an image; 
   *  The warped image defines the geometry and will contain the resampled image contents. 
   *  \param sourceImage The source image will be resampeled into the warped image space according to the transformation parameters
   *  \param warpedImage The warped image defines the geometry of the internally calculated DVF and will also hold the actually warped image.
   */
  void TransformImage( nifti_image* sourceImage, nifti_image* warpedImage );

  /** Transform an image using the "push-interpolation". This function was originally implemented
   *  specifically for reg_resp, hence the functionality will be provided by protected function of 
   *  this class.
   *  \param sourceImage Image that will be push-transformed into the warped image
   *  \param sourceWeightsImage Weights associated with the source image. This considers 
   *                            the image acquisition process where the acquired image 
   *                            may not have contributed equally to the voxel in the source 
   *                            image. 
   *  \param warpedImage Image into which the intensities will be pushed. Note that 
   *                     intensities will be accumulated.
   *  \param warpedWeightsImage Warped (and accumulated) weights of the source image.
   */
  void TransformImageAdjoint( nifti_image* sourceImage, nifti_image* sourceWeightsImage,
    nifti_image* warpedImage, nifti_image* warpedWeightsImage );

  /** Generate a deep copy of the transformation and return it.
   */
  virtual std::shared_ptr<Transformation> DeepCopy() = 0;

  /** Print the current transformation parameters to the command line.
   */
  virtual void DisplayTransformationParameters() = 0;

  /** Calculate the sum of the penalty weights
   */
  virtual PrecisionType GetSumOfPenaltyWeights() = 0;

  /** Calculate the gradient of an image with respect to the deformation vector field.
   *  This calculates \f$ \frac{\partial \mathbf{I}_{T_t}}{\partial \textbf{DVF}_t} \f$ and is
   *  not dependent on any parametrisation of a derived transformation. 
   *  \param sourceImage Nifti-image structure with the image from which the derivative will be calculated.
   *  \param outWarpedGradientImage Nifti image structure wich will be filled with the calculated gradient. 
                                    Has to be allocated outside of this function.
   */
  virtual void GetImageGradientWRTDVF(nifti_image* sourceImage, nifti_image* outWarpedGradientImage);

  /** Function to reorientate vectors in 5D nifti image.
   *  Vectors which are stored along the 5th (u) dimension will be reoriented
   *  \param vectorFieldToReorientate Vector image to be reorientated
   *  \param reorientationMatrix 4 x 4 Matrix used to reorientate vectors
   */
  void ReorientateVectorImage(nifti_image* vectorFieldToReorientate, mat44 reorientationMatrix);


  /** Function to determine the maximum length of a given set of parameters. 
   *   This length depends on how the transformation parameters correspond internally to each other, i.e. how the x, y, and z
   *   component are stored.
   * 
   */
  virtual PrecisionType GetMaxTransformationParameterLength( PrecisionType* parametersIn ) = 0;


  /** Function that returns the current transformation a vecotr of images. These are not the DVFs, but images of the parameters instead. 
   *  Hence it can only be implemented by any derived classes. 
   */
  virtual std::vector<nifti_image*> GetTransformationAsImage() = 0;


  /** Get the padding value.
   */
  inline PrecisionType GetPaddingValue() { return this->warpedPaddingValue; };

protected:

  /** Function to check if the DVF needs to be updated according to the target image
   *  \param targetImageIn The geometry if the DVF is checked against the target image
   */
  bool CheckDVFImageUpdateRequired( nifti_image* targetImageIn );


  /** Function that checks the basic geometry of two images. Is used by CheckDVFImageUpdateRequired() and 
   *  checks the image size (nx, ny, and nz) as well as the qform/sform matrix depending on the qform/sform 
   *  code. 
   *  Note: Differences in higher dimensions (nu, nv, nt) are allowed to facilitd comparison between target 
   *  image and DVF. 
   */
  bool CheckImageGeometryEquality( nifti_image* img1, nifti_image* img2 );


  /** Helper function to delete the DVF image
   */
  virtual void ClearDeformationVectorFieldImage();

  // Functions brought forward from reg-resp
  void CubicSplineTransformImage3D(nifti_image* sourceImage, nifti_image* warpedImage);
  void CubicSplineTransformImage2D( nifti_image* sourceImage, nifti_image* warpedImage );
  void NearestNeighbourTransformImage3D( nifti_image* sourceImage, nifti_image* warpedImage );
  void NearestNeighbourTransformImage2D( nifti_image* sourceImage, nifti_image* warpedImage );
  void TrilinearTransformImage( nifti_image* sourceImage, nifti_image* warpedImage );
  void BilinearTransformImage( nifti_image* sourceImage, nifti_image* warpedImage );

  void CubicSplineTransformImageAdjoint3D( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );
  void CubicSplineTransformImageAdjoint2D( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );
  void NearestNeighbourTransformImageAdjoint3D( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );
  void NearestNeighbourTransformImageAdjoint2D( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );
  void TrilinearTransformImageAdjoint( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );
  void BilinearTransformImageAdjoint( nifti_image * sourceImage, nifti_image * sourceWeightsImage, nifti_image * warpedImage, nifti_image * warpedWeightsImage );

  void TrilinearImageGradient( nifti_image* sourceImage, nifti_image* resultGradientImage );
  void BilinearImageGradient( nifti_image* sourceImage, nifti_image* resultGradientImage );
  void CubicSplineImageGradient3D( nifti_image* sourceImage, nifti_image* resultGradientImage );
  void CubicSplineImageGradient2D( nifti_image* sourceImage, nifti_image* resultGradientImage );

  bool dvfImageUpdateRequired;                 ///< Indicates if the DVF requires re-calculation when it is required
  unsigned int numberOfParameters;             ///< The total number of parameters needed to describe the transformation
  int interpolation;                           ///< Type of interpolation used \todo Change into enumeration.
  PrecisionType warpedPaddingValue;            ///< Padding value for the warped image
  nifti_image* deformationVectorFieldImage;    ///< The deformation vector field image
};

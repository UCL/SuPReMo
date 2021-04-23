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

//----------
// includes 
//----------
#include <memory>
#include <vector>
#include "nifti1_io.h"
#include "SupremoUtils.h"



//----------------------
// forward declarations
//----------------------
class CorrespondenceModel;
class ImageAcquisition;




/** Baseclass defining the interface and some default implementations for the motion-compensated image reconstruction.  
 */
class MoCoRecon
{
public:
  typedef float PrecisionType;  
  typedef std::vector< std::vector<PrecisionType> > SurrogateSignalType;

  /** Constructor
   */
  MoCoRecon();

  /** Destructor
   */
  virtual ~MoCoRecon();

  /** Set the correspondence model. 
   *  \param correspondenceModelIn The input correspondence model
   */
  void SetCorrespondenceModel( const std::shared_ptr<CorrespondenceModel>& correspondenceModelIn );
  
  /** Set the image acquisition object. 
   *  \param imageAcquisitionIn The object providing the image-acquisition functionality.
   */
  void SetImageAcquisition( const std::shared_ptr<ImageAcquisition>& imageAcquisitionIn );

  /** Set the surrogate signal. Has to be given, such that the correspondence model can be used to 
   *  generate the required transformations.
   */
  void SetSurrogateSignals( const SurrogateSignalType & surrSignalsIn );
  
  /** Define the dynamic images
   *  \param Vector holding the pointer to the dynamic images. 
   */
  void SetDynamicImages( const std::vector<nifti_image*>& dynamicImagesIn );
  
  /** Define the reconstruction geometry. 
   *  \param reconstructionGeometryImageIn Pointer to the image that defines the geometry of the reconstructed image. 
   */
  void SetReconstructionGeometryImage( nifti_image* reconstructionGeometryImageIn );
  
  /** Abstract base class to actually calculate the motion-compensated reconstruction. Has to be implemented by a specific method. 
   */
  virtual void Update() = 0;

  /** Get a pointer to the reconstructed image. Note, the destructor of this class will clean up the image. 
   */
  nifti_image* GetReconstructedImage();


  /** Returns a boolean indicating if a repeated call of the Update() method changes the result - even if 
   *  the CorrespondenceModel, or the dynamic images did not change. This behaviour has to be determined by a 
   *  derived class by setting repeatedUpdateChangesMCRImage accordingly.
   */
  bool GetRepeatedUpdateChangesMCR() { return this->repeatedUpdateChangesMCRImage; };

  /** Copy the reconstructed image to the destination image. 
   *  \param destinationImage Nifti-image into which the reconstructed image will be written. Has to have the same size as the reconstruction geometry image.
   */
  void CopyReconstructedImageContentsToImage( nifti_image* destinationImage );

protected:
  /** Allocate memory to hold the reconstructed image
   */
  void AllocateReconstructedImage();

  /** Clear the memory reserved for the reconstructed image and set the corresponding pointer to null.
   */
  void ClearReconstructedImage();

  std::shared_ptr<CorrespondenceModel> correspondenceModel;  ///< Pointer to the correspondence model 
  std::shared_ptr<ImageAcquisition> imageAcquisition;        ///< Pointer to the image acquisition opbject 
  SurrogateSignalType surrogateSignals;                      ///< Vector holding the surrogate signals
  std::vector<nifti_image*> dynamicImages;                   ///< Vector holding the dynamic images
  size_t numberOfDynamicImages;                              ///< The number of the dynamic images - will be derived from the vector dynamicImages
  nifti_image* reconstructionGeometryImage;                  ///< Pointer to the image that defines the reconstruction geometry, think of it as a reference image
  nifti_image* reconstructedImage;                           ///< Pointer to the reconstructed image
  PrecisionType paddingValue;                                ///< The padding value used if no data available to reconstruct an image intensity value
  bool repeatedUpdateChangesMCRImage;                        ///< Indicates wether a repeated call of the update function cahgnes the MCR result. 
};

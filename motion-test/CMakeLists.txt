#  ====================================================================================================   
#                                                                                                         
#    SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
#             An implementation of the generalised motion modelling and image registration framework      
#                                                                                                         
#    Copyright (c) University College London (UCL). All rights reserved.                                  
#                                                                                                         
#    This software is distributed WITHOUT ANY WARRANTY; without even                                      
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
#    PURPOSE.                                                                                             
#                                                                                                         
#    See LICENSE.txt in the top level directory for details.                                              
#                                                                                                         
#  ====================================================================================================   




# Set the data folder 
set( DATAFOLDER ${CMAKE_CURRENT_SOURCE_DIR}/test-data )

# Set the include directory
include_directories(${PROJECT_SOURCE_DIR}/motion-lib ${PROJECT_SOURCE_DIR}/motion-ext)

# Generate the executable and link 
add_executable(testImagePyramidGeneration testImagePyramidGeneration.cpp)
target_link_libraries(testImagePyramidGeneration ${NiftyReg_LIBRARIES} )

add_executable(testSSDImageSimilarity testSSDImageSimilarity.cpp)
target_link_libraries(testSSDImageSimilarity ImageSimilarity ${NiftyReg_LIBRARIES} )

add_executable(testBSplineTransformation testBSplineTransformation.cpp)
target_link_libraries(testBSplineTransformation Transformation ${NiftyReg_LIBRARIES} )

add_executable(testSlidingBSplineTransformation testSlidingBSplineTransformation.cpp)
target_link_libraries(testSlidingBSplineTransformation Transformation CommandLineParser ${NiftyReg_LIBRARIES} )

add_executable(testAdjointTransformation testAdjointTransformation.cpp)
target_link_libraries(testAdjointTransformation Transformation CommandLineParser ${NiftyReg_LIBRARIES} )

add_executable(testObjectiveFunction testObjectiveFunction.cpp)
target_link_libraries(testObjectiveFunction Transformation ImageSimilarity CorrespondenceModel CommandLineParser ObjectiveFunction ImageAcquisition ${NiftyReg_LIBRARIES} )

add_executable(testWeightedAverageMoCoReco testWeightedAverageMoCoReco.cpp)
target_link_libraries(testWeightedAverageMoCoReco Transformation CorrespondenceModel CommandLineParser MoCoRecon ImageAcquisition ${NiftyReg_LIBRARIES} )

add_executable(testSuperResIBPMoCoRecoRestart testSuperResIBPMoCoRecoRestart.cpp)
target_link_libraries(testSuperResIBPMoCoRecoRestart Transformation CorrespondenceModel CommandLineParser MoCoRecon ImageAcquisition ${SUPREMO_EXT_LIBS} ${NiftyReg_LIBRARIES} )

add_executable(testSuperResIBPMoCoRecoUpdate testSuperResIBPMoCoRecoUpdate.cpp)
target_link_libraries(testSuperResIBPMoCoRecoUpdate Transformation CorrespondenceModel CommandLineParser MoCoRecon ImageAcquisition ${SUPREMO_EXT_LIBS} ${NiftyReg_LIBRARIES} )

add_executable(testLowResolutionImageAcquisition testLowResolutionImageAcquisition.cpp)
target_link_libraries(testLowResolutionImageAcquisition CommandLineParser ImageAcquisition ${NiftyReg_LIBRARIES})

add_executable(testLowResolutionImageAcquisitionAllocMinSizeImg testLowResolutionImageAcquisitionAllocMinSizeImg.cpp)
target_link_libraries(testLowResolutionImageAcquisitionAllocMinSizeImg  CommandLineParser ImageAcquisition ${NiftyReg_LIBRARIES})

add_executable(testLowResolutionImageAcquisitionAdjoint testLowResolutionImageAcquisitionAdjoint.cpp)
target_link_libraries( testLowResolutionImageAcquisitionAdjoint CommandLineParser ImageAcquisition ${NiftyReg_LIBRARIES} )


# Add the test to be performed
# Image pyramids
add_test( ImagePyramidTest_2D_L0 testImagePyramidGeneration ${DATAFOLDER}/2d_slice_dyn22.nii.gz 3 3 0 ${DATAFOLDER}/expectedOutcomes/2d_slice_dyn22_pyramid_L0.nii.gz)
add_test( ImagePyramidTest_2D_L1 testImagePyramidGeneration ${DATAFOLDER}/2d_slice_dyn22.nii.gz 3 3 1 ${DATAFOLDER}/expectedOutcomes/2d_slice_dyn22_pyramid_L1.nii.gz)
add_test( ImagePyramidTest_2D_L2 testImagePyramidGeneration ${DATAFOLDER}/2d_slice_dyn22.nii.gz 3 3 2 ${DATAFOLDER}/expectedOutcomes/2d_slice_dyn22_pyramid_L2.nii.gz)

add_test( ImagePyramidTest_3D_L0 testImagePyramidGeneration ${DATAFOLDER}/3d_vol_static.nii.gz 3 3 0 ${DATAFOLDER}/expectedOutcomes/3d_vol_pyramid_L0.nii.gz)
add_test( ImagePyramidTest_3D_L1 testImagePyramidGeneration ${DATAFOLDER}/3d_vol_static.nii.gz 3 3 1 ${DATAFOLDER}/expectedOutcomes/3d_vol_pyramid_L1.nii.gz)
add_test( ImagePyramidTest_3D_L2 testImagePyramidGeneration ${DATAFOLDER}/3d_vol_static.nii.gz 3 3 2 ${DATAFOLDER}/expectedOutcomes/3d_vol_pyramid_L2.nii.gz)

# Image similarity
add_test( ImageSimilarityMeasure_2D testSSDImageSimilarity ${DATAFOLDER}/test_simil_ref_image_22.nii.gz ${DATAFOLDER}/test_simil_simDyn_image_22.nii.gz ${DATAFOLDER}/expectedOutcomes/ssd_measure_22.txt )

# BSpline transformation
add_test( BSplineTransformation_2D testBSplineTransformation ${DATAFOLDER}/test_warpWithCPG_cpg_22.nii.gz ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz ${DATAFOLDER}/test_warpWithCPG_a_reference_22_f.nii.gz            ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz ${DATAFOLDER}/expectedOutcomes/test_warpWithCPG_a_dvf_22.nii.gz ${DATAFOLDER}/expectedOutcomes/test_warpWithCPG_a_warped_22.nii.gz )
add_test( BSplineTransformation_3D testBSplineTransformation ${DATAFOLDER}/test_warpWithCPG_cpg_22.nii.gz ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz ${DATAFOLDER}/expectedOutcomes/test_warpWithCPG_b_dvf_22.nii.gz ${DATAFOLDER}/expectedOutcomes/test_warpWithCPG_b_warped_22.nii.gz )

# Adjoint transformation (B-spline based)
add_test( NAME AdjointTransformationTest COMMAND testAdjointTransformation
			-cpg    ${DATAFOLDER}/test_warpWithCPG_cpg_22.nii.gz 
			-cpgRef ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22.nii.gz 
			-imgX   ${DATAFOLDER}/test_warpWithCPG_ab_floating_b_reference22_toZero.nii.gz 
			-imgY   ${DATAFOLDER}/expectedOutcomes/test_warpWithCPG_b_warped_22_toZero.nii.gz)

# Sliding B-spline transformation 
add_test( NAME SlidingBSplineTransformation_3D_DVF COMMAND testSlidingBSplineTransformation
      -cpg1 ${DATAFOLDER}/xcat/sliding/cpp_mat1to6Sliding_region1.nii.gz
      -cpg2 ${DATAFOLDER}/xcat/sliding/cpp_mat1to6Sliding_region2.nii.gz
      -distMap ${DATAFOLDER}/xcat/sliding/distMap1.nii.gz
      -source ${DATAFOLDER}/xcat/mat_1.nii.gz
      -target ${DATAFOLDER}/xcat/mat_6.nii.gz
      -expectedDVF ${DATAFOLDER}/expectedOutcomes/sliding/dvf_combined.nii.gz
      -expectedTransSource ${DATAFOLDER}/expectedOutcomes/sliding/warped.nii.gz
 )

add_test( NAME SlidingBSplineTransformation_3D_GOCT_GRAD COMMAND testSlidingBSplineTransformation
      -cpg1 ${DATAFOLDER}/xcat/sliding/cpp_mat1to6Sliding_region1.nii.gz
      -cpg2 ${DATAFOLDER}/xcat/sliding/cpp_mat1to6Sliding_region2.nii.gz
      -distMap ${DATAFOLDER}/xcat/sliding/distMap1.nii.gz
      -source ${DATAFOLDER}/xcat/mat_1.nii.gz
      -target ${DATAFOLDER}/xcat/mat_6.nii.gz
      -expectedGOCTGrad1 ${DATAFOLDER}/expectedOutcomes/sliding/goctOnly_w0.01_region1.nii.gz
      -expectedGOCTGrad2 ${DATAFOLDER}/expectedOutcomes/sliding/goctOnly_w0.01_region2.nii.gz
 )

# Objective function
add_test( NAME ObjectiveFunctionTest COMMAND testObjectiveFunction 
	        -refState ${DATAFOLDER}/xcat/mat_1.nii.gz 
	        -dynamic 10  ${DATAFOLDER}/xcat/volumes_XCAT_fileList.txt
	        -surr 2 ${DATAFOLDER}/xcat/surrogate_XCAT_ApSi_Zsc.txt
	        -valFile ${DATAFOLDER}/expectedOutcomes/objectiveFunctionValue_at0.txt
	        -gradFiles ${DATAFOLDER}/expectedOutcomes/test_rcm_grad_lev_0_iter_0_surr_0.nii.gz,${DATAFOLDER}/expectedOutcomes/test_rcm_grad_lev_0_iter_0_surr_1.nii.gz
	        WORKING_DIRECTORY ${DATAFOLDER}  )

# Weighted average motion-compensated image reconstruction
add_test( NAME WeightedAverageMoCoRecoTest COMMAND  testWeightedAverageMoCoReco 
	        -refState ${DATAFOLDER}/xcat/mat_1.nii.gz 
	        -dynamic 10  ${DATAFOLDER}/xcat/volumes_XCAT_fileList.txt
	        -surr 2 ${DATAFOLDER}/xcat/surrogate_XCAT_ApSi_Zsc.txt
	 		-rcmIn ${DATAFOLDER}/respRCM_fitIt20_switchIt10_l3_sx10_be01.nii.gz 
	 		-mcrToCompare ${DATAFOLDER}/expectedOutcomes/respMCR_fitIt20_switchIt10_l3_sx10_be01.nii.gz 
			-ln 3
	 		WORKING_DIRECTORY ${DATAFOLDER} )

# Super resolution IBP motion-compensated image reconstruction (rcm at start)
add_test( NAME SuperResolutionIBPRestartMoCoReco_initial_Test COMMAND  testSuperResIBPMoCoRecoRestart 
              -refState ${DATAFOLDER}/xcat/MRsim/target_2x2x2.nii.gz 
              -dynamic 355  ${DATAFOLDER}/xcat/MRsim/modelSlices_tumourA_sagAxi_R1.txt
              -surr 2 ${DATAFOLDER}/xcat/MRsim/surrSkinDiaDetected_tumourA_sagAxi_R1.txt
              -mcrIt 2
              -ln 3
              -rcmIn ${DATAFOLDER}/xcat/MRsim/RCM_lev_0_initial_type2_itNum2.nii.gz
              -mcrToCompare ${DATAFOLDER}/expectedOutcomes/MCR_IBP_lev_0_initial_type2_itNum2.nii.gz
              WORKING_DIRECTORY ${DATAFOLDER} )


# Super resolution IBP motion-compensated image reconstruction (rcm after first iteration)
add_test( NAME SuperResolutionIBPRestartMoCoReco_firstIt_Test COMMAND  testSuperResIBPMoCoRecoRestart 
              -refState ${DATAFOLDER}/xcat/MRsim/target_2x2x2.nii.gz 
              -dynamic 355  ${DATAFOLDER}/xcat/MRsim/modelSlices_tumourA_sagAxi_R1.txt
              -surr 2 ${DATAFOLDER}/xcat/MRsim/surrSkinDiaDetected_tumourA_sagAxi_R1.txt
              -mcrIt 2
              -ln 3
              -rcmIn ${DATAFOLDER}/xcat/MRsim/RCM_lev_0_iter_0_type2_itNum2.nii.gz
              -mcrToCompare ${DATAFOLDER}/expectedOutcomes/MCR_IBP_lev_0_iter_0_type2_itNum2.nii.gz
              WORKING_DIRECTORY ${DATAFOLDER} )

# Super resolution IBP motion-compensated image reconstruction (rcm after first iteration, updating reconstruction)
add_test( NAME SuperResolutionIBPUpdateMoCoReco_firstIt_Test COMMAND  testSuperResIBPMoCoRecoUpdate 
              -refState ${DATAFOLDER}/xcat/MRsim/target_2x2x2.nii.gz 
              -dynamic 355  ${DATAFOLDER}/xcat/MRsim/modelSlices_tumourA_sagAxi_R1.txt
              -surr 2 ${DATAFOLDER}/xcat/MRsim/surrSkinDiaDetected_tumourA_sagAxi_R1.txt
              -mcrIt 1
              -ln 3
              -rcmIn ${DATAFOLDER}/xcat/MRsim/RCM_lev_0_iter_0_type3_itNum1.nii.gz
              -mcrToCompare ${DATAFOLDER}/expectedOutcomes/MCR_IBP_lev_0_iter_0_type3_itNum1.nii.gz
              WORKING_DIRECTORY ${DATAFOLDER} )

# Low-resolution image acquisition
add_test( NAME LowResolutionImageAcquisitionSimulationSAG COMMAND testLowResolutionImageAcquisition 
            -dyn ${DATAFOLDER}/lowRes_curReference_image_11_ln_0_of_2.nii.gz 
            -fullImg ${DATAFOLDER}/lowRes_warped_image_11_ln_0_of_2.nii.gz 
            -resSimLow ${DATAFOLDER}/expectedOutcomes/lowRes_simDyn_image_11_ln_0_of_2.nii.gz 
            WORKING_DIRECTORY ${DATAFOLDER} )

add_test( NAME LowResolutionImageAcquisitionSimulationAXI COMMAND testLowResolutionImageAcquisition 
            -dyn ${DATAFOLDER}/lowRes_curReference_image_197_ln_0_of_2.nii.gz
            -fullImg ${DATAFOLDER}/lowRes_warped_image_197_ln_0_of_2.nii.gz
            -resSimLow ${DATAFOLDER}/expectedOutcomes/lowRes_simDyn_image_197_ln_0_of_2.nii.gz
            WORKING_DIRECTORY ${DATAFOLDER} )

add_test( NAME LowResolutionImageAcquisitionSimulationMinSizeImgAllocSAG COMMAND testLowResolutionImageAcquisitionAllocMinSizeImg 
            -dyn ${DATAFOLDER}/lowRes_curReference_image_11_ln_0_of_2.nii.gz 
            -fullImg ${DATAFOLDER}/lowRes_curFloating_image_11_ln_0_of_2.nii.gz 
            -resMinSizeImg ${DATAFOLDER}/lowRes_warped_image_11_ln_0_of_2.nii.gz 
            WORKING_DIRECTORY ${DATAFOLDER} )

add_test( NAME LowResolutionImageAcquisitionSimulationMinSizeImgAllocAXI COMMAND testLowResolutionImageAcquisitionAllocMinSizeImg 
            -dyn ${DATAFOLDER}/lowRes_curReference_image_197_ln_0_of_2.nii.gz 
            -fullImg ${DATAFOLDER}/lowRes_curFloating_image_197_ln_0_of_2.nii.gz 
            -resMinSizeImg ${DATAFOLDER}/lowRes_warped_image_197_ln_0_of_2.nii.gz 
            WORKING_DIRECTORY ${DATAFOLDER} )

add_test( NAME testLowResolutionImageAcquisitionAdjointSAG COMMAND testLowResolutionImageAcquisitionAdjoint
            -dyn ${DATAFOLDER}/lowRes_curReference_image_11_ln_0_of_2.nii.gz 
            -fullImg ${DATAFOLDER}/lowRes_curFloating_image_11_ln_0_of_2.nii.gz 
            -resAdjImg ${DATAFOLDER}/expectedOutcomes/lowRes_dynamicAfterAdjoint_image_11_ln_0_of_2.nii.gz
            -resAdjWeightImg ${DATAFOLDER}/expectedOutcomes/lowRes_weightsInDynAfterAdjImage_image_11_ln_0_of_2.nii.gz
            WORKING_DIRECTORY ${DATAFOLDER} )

add_test( NAME testLowResolutionImageAcquisitionAdjointAXI COMMAND testLowResolutionImageAcquisitionAdjoint
            -dyn ${DATAFOLDER}/lowRes_curReference_image_197_ln_0_of_2.nii.gz 
            -fullImg ${DATAFOLDER}/lowRes_curFloating_image_197_ln_0_of_2.nii.gz 
            -resAdjImg ${DATAFOLDER}/expectedOutcomes/lowRes_dynamicAfterAdjoint_image_197_ln_0_of_2.nii.gz
            -resAdjWeightImg ${DATAFOLDER}/expectedOutcomes/lowRes_weightsInDynAfterAdjImage_image_197_ln_0_of_2.nii.gz
            WORKING_DIRECTORY ${DATAFOLDER} )


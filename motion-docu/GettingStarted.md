
# Getting started with SuPReMo {#gettingstarted}


## Building from source

## Dependencies

### CMake

The cross platform make-file generator is used to manage the build on various platforms. 
CMake can be downloaded [here](https://cmake.org). Please install a reasonably recent version 
on your system, a minimum of version 3.3 is required by the CMakeLists.txt. 

### Nifty-Reg

During the build step SuPReMo links to the libraries of *nifty-reg* which is available 
[here](https://github.com/KCL-BMEIS/niftyreg) and the build will be described below.


## Building SuPReMo on a Linux system

The steps below describe how to build SuPReMo (and nifty-reg) on a Linux system. We perform an out-of source 
build for both, nifty-reg and SuPReMo. Choose a location where to download and build the software packages, 
we will name this `WORK_DIR` in the instructions below.

### Building nifty-reg

To build the required *nifty-reg* libraries in the `WORK_DIR` download the source code for instacne using git, 
create a build and an install directory and run `ccmake` or `cmake-gui` in the build directory
to configure the build.  

```bash
 cd WORK_DIR
 git clone --single-branch --branch 28-sliding-regions https://github.com/KCL-BMEIS/niftyreg
 mkdir niftyreg_build
 mkdir niftyreg_install
 cd niftyreg_build
 ccmake ../niftyreg
```

**Note:** Using the branch `28-sliding-regions` fixes some open-mp errors that occured with gcc version 9, however, gcc version 7 works well with 
the master branch. When using an older compiler you can chose to replace the above git-clone with the simpler one below. The branch sliding-regions is
not a pre-requisite to utilise the sliding implementation of SuPReMo:
```bash
 git clone https://github.com/KCL-BMEIS/niftyreg

```

Configure the nifty-reg build (by pressing `c`) and set the following options within ccmake using the arrow keys
to move up and down the options, select one by pressing `enter` and then type the required option value and confirm 
with `enter`. When all parameters are set, press `c` and then `g` to generate the `Makefile`. The cmake-gui is potentially more convenient to 
use but the result will be the same. 

- `BUILD_ALL_DEP        = ON`
- `BUILD_SHARED_LIBS    = OFF`
- `BUILD_TESTING        = OFF`
- `CMAKE_BUILD_TYPE     = Release`
- `CMAKE_INSTALL_PREFIX = WORK_DIR/niftyreg_install` 
- `USE_CUDA             = OFF`
- `USE_OPENCL           = OFF`
- `USE_OPENMP           = ON`
- `USE_SSE              = ON`

Still within the build-directory start the build and install the binaries to the directory defined in the 
configuration step above (`WORK_DIR/niftyreg_install`)
```bash
make
make install

```
to start the build process. On a multi-core system you may speed up the build by allowing parallel build steps where 
possible by usint the option `-j4`, for four parallel builds for instance. 

To check that the installation was successful, you can get a help message from nifty-reg:
```bash
cd WORK_DIR/niftyreg_install/bin/runSupremo
reg_f3d -h

```
If a detailed usage message shows up on the command line, you are good to proceed. 

### Building SuPReMo {#buildingsupremo}

The build steps for SuPReMo are very similar to the ones required for nifty-reg, namely, obtaining the source code, configuring the build, building and installing.


```bash
 cd WORK_DIR
 git clone https://github.com/UCL/SuPReMo
 mkdir SuPReMo_build
 mkdir SuPReMo_install
 cd SuPReMo_build
 ccmake ../SuPReMo

```

Configure SuPReMo with the following optoins. Note however, that since nifty-reg is a dependency of SuPReMo, on the first connfiguration attempt 
an error will be raised, simply because the install path of nifty-reg is not known to cmake at this stage. Providing it during the configuration
solves this. For this set the variable NiftyReg_DIR to `WORK_DIR/niftyreg_install` and type `c`. The variables `NiftyReg_INCLUDE_DIR` and 
`NiftyReg_TOOLS_LIBRARY` will be filled in automatically.
- `BUILD_TESTING        = ON/OFF` depending on if you want to run the tests on your local system
- `CMAKE_BUILD_TYPE     = Release`
- `CMAKE_INSTALL_PREFIX = WORK_DIR/SuPReMo_install/`
- `NiftyReg_DIR         = WORK_DIR/niftyreg_install/`
- `USE_OPENMP           = ON`

After the configuration a message will be displayed saying *niftyreg found*. This must be acknowledged wtith `e`. 

After the configuraiton, SuPReMo can be build and installed:
```bash
make
make install

```

If successful you will have the installation in `WORK_DIR/SuPReMo_install` and the main executable in the folder `bin`. You display SuPReMo's
usage message on the command line by runing
```bash
WORK_DIR/SuPReMo_install/bin/runSupremo -h

```

### Running the tests included in SuPReMo

If you want to check your local SuPReMo installation more thoroughly, you can run the tests using ctest. Navigate to the build-directory 
and call ctest as follows:
```bash
ctest --output-on-failure -C Release

```
Which will generate an output similar to this:
```
Test project WORK_DIR/SuPReMo_build
      Start  1: ImagePyramidTest_2D_L0
 1/23 Test  #1: ImagePyramidTest_2D_L0 ......................................   Passed    0.11 sec
      Start  2: ImagePyramidTest_2D_L1
 2/23 Test  #2: ImagePyramidTest_2D_L1 ......................................   Passed    0.10 sec
      Start  3: ImagePyramidTest_2D_L2
 3/23 Test  #3: ImagePyramidTest_2D_L2 ......................................   Passed    0.09 sec
      Start  4: ImagePyramidTest_3D_L0
 4/23 Test  #4: ImagePyramidTest_3D_L0 ......................................   Passed    0.15 sec
      Start  5: ImagePyramidTest_3D_L1
 5/23 Test  #5: ImagePyramidTest_3D_L1 ......................................   Passed    0.15 sec
      Start  6: ImagePyramidTest_3D_L2
 6/23 Test  #6: ImagePyramidTest_3D_L2 ......................................   Passed    0.18 sec
      Start  7: ImageSimilarityMeasure_2D
 7/23 Test  #7: ImageSimilarityMeasure_2D ...................................   Passed    0.08 sec
      Start  8: BSplineTransformation_2D
 8/23 Test  #8: BSplineTransformation_2D ....................................   Passed    0.26 sec
      Start  9: BSplineTransformation_3D
 9/23 Test  #9: BSplineTransformation_3D ....................................   Passed    0.49 sec
      Start 10: AdjointTransformationTest
10/23 Test #10: AdjointTransformationTest ...................................   Passed    0.32 sec
      Start 11: SlidingBSplineTransformation_3D_DVF
11/23 Test #11: SlidingBSplineTransformation_3D_DVF .........................   Passed    2.48 sec
      Start 12: SlidingBSplineTransformation_3D_GOCT_GRAD
12/23 Test #12: SlidingBSplineTransformation_3D_GOCT_GRAD ...................   Passed    1.64 sec
      Start 13: ObjectiveFunctionTest
13/23 Test #13: ObjectiveFunctionTest .......................................   Passed    9.58 sec
      Start 14: WeightedAverageMoCoRecoTest
14/23 Test #14: WeightedAverageMoCoRecoTest .................................   Passed    7.29 sec
      Start 15: SuperResolutionIBPRestartMoCoReco_initial_Test
15/23 Test #15: SuperResolutionIBPRestartMoCoReco_initial_Test ..............   Passed    8.32 sec
      Start 16: SuperResolutionIBPRestartMoCoReco_firstIt_Test
16/23 Test #16: SuperResolutionIBPRestartMoCoReco_firstIt_Test ..............   Passed    7.28 sec
      Start 17: SuperResolutionIBPUpdateMoCoReco_firstIt_Test
17/23 Test #17: SuperResolutionIBPUpdateMoCoReco_firstIt_Test ...............   Passed    8.09 sec
      Start 18: LowResolutionImageAcquisitionSimulationSAG
18/23 Test #18: LowResolutionImageAcquisitionSimulationSAG ..................   Passed    0.07 sec
      Start 19: LowResolutionImageAcquisitionSimulationAXI
19/23 Test #19: LowResolutionImageAcquisitionSimulationAXI ..................   Passed    0.08 sec
      Start 20: LowResolutionImageAcquisitionSimulationMinSizeImgAllocSAG
20/23 Test #20: LowResolutionImageAcquisitionSimulationMinSizeImgAllocSAG ...   Passed    0.16 sec
      Start 21: LowResolutionImageAcquisitionSimulationMinSizeImgAllocAXI
21/23 Test #21: LowResolutionImageAcquisitionSimulationMinSizeImgAllocAXI ...   Passed    0.15 sec
      Start 22: testLowResolutionImageAcquisitionAdjointSAG
22/23 Test #22: testLowResolutionImageAcquisitionAdjointSAG .................   Passed    0.19 sec
      Start 23: testLowResolutionImageAcquisitionAdjointAXI
23/23 Test #23: testLowResolutionImageAcquisitionAdjointAXI .................   Passed    0.20 sec

100% tests passed, 0 tests failed out of 23

Total Test time (real) =  47.75 sec

```



\tableofcontents

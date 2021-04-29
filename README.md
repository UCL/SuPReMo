[![Build Status](https://travis-ci.com/UCL/SuPReMo.svg?token=ryYEyvSMfhuCX4W6bRRD&branch=main)](https://travis-ci.com/UCL/SuPReMo)

# SuPReMo

Software implementing the unified motion modelling and image registation 
framework as described by McClelland et al. 

# Transition from nifty-reg resp

**Note**: The software is still in its early stages and may show some unusual 
or unwanted behaviour. Please report this by either emailing me or creating an 
issue. Thanks!


Parameters: 

| reg_resp       | runSupremo    | Pattern       | Description                                                                         |
|----------------|---------------|---------------|-------------------------------------------------------------------------------------|
| -static        | -refState     | [fName]       | Reference state/static image                                                        |
| -dynamic       | -dynamic      | [int] [fName] | File containing list of dynamic images and the total number of dynamic images       |
| -surr          | -surr         | [int] [fName] | File containting the surrogate signal and the number of signals used                |
| -dType         | -dType        | [int]         | Dynamic image data type [0] = full res images, [1] = low res images                 |
| -defSpace      | -defSpace     | [fName]       | File name of image that defines the space of the deformed iamges                    |
| -mcrType       | -mcrType      | [int]         | Motion compensation [0] none, [1] average weighting, [2] SR, restart, [3] SR update |
| -maxMCRIt      | -maxMCRIt     | [int]         | Number of iterations used for iterative MCR method                                  |
| -rcm           | -outRCM       | [fName]       | Path to save the respratory correspondence model (RCM) to                           |
| -saveDyn       | -outSimDyn    | [fName]       | Path to save the simulated dynamic images to                                        |  
| -saveMCR       | -outMCR       | [fName]       | Path to save the motion-compensated reconstruction images to                        |  
| -saveInterMCRs | -outInterMCR  | [fName]       | Path to save the motion-compensated reconstruction images to                        |  
|                | -outInterGrad | [fName]       | Path to save the intermediate gradient to (experimental)                            |  
| -sx            | -sx           | [float]       | Final grid spacing along x axis (in mm if positive, in voxels if negative)          |
| -sy            | -sy           | [float]       | Final grid spacing along y axis (in mm if positive, in voxels if negative)          |  
| -sz            | -sz           | [float]       | Final grid spacing along z axis (in mm if positive, in voxels if negative)          |  
| -be            | -be           | [float]       | Bending energy constraint weight                                                    |  
|                | -le           | [float]       | Linear energy constraint weight (experimental)                                      |  
| -maxSwitchIt   | -maxSwitchIt  | [int]         | Maximum number of switches between fitting and reconstruction                       |
| -ln            | -ln           | [int]         | Number of pypramid levels generated                                                 |  
| -lp            | -lp           | [int]         | Number of pypramid levels used                                                      |  
| -maxFitIt      | -maxFitIt     | [int]         | Number of iterations to fit the motion model                                        |  

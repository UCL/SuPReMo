# Background {#background}

The theoretical background is explained in the paper by [Jamie McClelland et al.](https://doi.org/10.1088/1361-6560/aa6070). 
Some more specifics which are relevant to the implementation of SuPReMo are provided on this page. 

## Nomenclature

| Symbol | Description |
|--------|-------------|
| \f$ \mathbf{I}_{t},\;t\in [1,\ldots, N_i] \f$ | motion images (or dynamic images, in image space) |
| \f$ t \f$| time point |
| \f$ N_i \f$ | total number of motion images | 
| \f$ \mathbf{I}_0 \f$ | reference state image | 
| \f$ \mathbf{M}_t = [m_{t,1}, \ldots, m_{t,N_m}]^\mathrm{T} \f$ |  motion parameters  |
| \f$ N_m \f$ | number of motion parameters |
| \f$ \mathbf{S}_t = [s_{t,1}, \ldots, s_{t,N_s} ]^\mathrm{T} \f$ | surrogate signal |
| \f$ N_s \f$ | number of surrogate signals per time point |
| \f$ \mathbf{R}=[\mathbf{R}_1,\ldots,\mathbf{R}_{N_r}]^\mathrm{T} \f$ | model parameters |
| \f$ N_r \f$ | number of model parameters per motion parameter |
| \f$ \mathbf{P}_{\negthickspace A_t} \f$ | simulated motion image (in acquisition space) |
| \f$ \mathbf{P}_{\negthickspace t}\f$ | motion image (in acquisition space) |
| \f$ N_{p} \f$ | number of voxels of the simulated dynamic image |

## Unifying image registration and respiratory motion modelling

The unified framework directly optimises the correspondence model parameters \f$\mathbf{R}\f$ on all motion images simultaneously. 

\anchor ObjectiveFunction
\f[
\mathcal{C}_t = (1-\lambda) \mathcal{S}( \mathbf{P}_{\negthickspace t}, A_t(\mathbf{T}(\mathbf{I}, \mathbf{M}_t)))
+ \lambda\mathcal{R}(\mathbf{M}_t)
\f]

Often the image similarity and regularisation are weighted with a parameter \f$\lambda \in ]0,1]\f$, such that a 
balance between data fidelity and smoothness of the solution can be controlled. 

In order to make use of gradient based optimisation methods, the partial derivative of the cost function \f$\mathcal{C}\f$ 
needs to be known with respect to those model parameters. 
\anchor DerivOfObjectiveFunctionWRTModelParams
\f[
\frac{\partial\mathcal{C}_t}{\partial \mathbf{R}} 
= 
\frac{\partial \mathbf{M}_t}{\partial\mathbf{R}}
\frac{\partial\mathcal{C}_t}{\partial\mathbf{M}_t}
\f]

### Correspondence model

The correspondence model relates the surrogate signal and the model parameters with the motion parameters, i.e.
\f[
\mathbf{M}_t = F(\mathbf{S}_t,\mathbf{R}).
\f]
For instance a 2nd order polynomial model can be expressed as
\anchor correspondenceModelQuadratic
\f[
\mathbf{M}_t = \mathbf{R}_1 s_t^2 +\mathbf{R}_2 s_t + \mathbf{R}_3
\f]
with the corresponding partial derivatives
\f{eqnarray*}{
\frac{\partial \mathbf{M}_t}{\partial \mathbf{R}_1} = s_t^2, \;\;
\frac{\partial \mathbf{M}_t}{\partial \mathbf{R}_2} = s_t, \;\;
\frac{\partial \mathbf{M}_t}{\partial \mathbf{R}_3} = 1.
\f}

## Algorithm

### Computation of objective function


The objective function value for a single motion-image time point is computed on the 
basis of \ref ObjectiveFunction "the objective function". The total objective function value is the sum over all time points
\anchor CostFuncSummedOverAllTimePoints
\f[
\mathcal{C}_\text{total} = \sum_{t=1}^{N_i} \mathcal{C}_t.
\f]


A procedural description of how the objective function value is calculated is presented in listing...
\todo Insert pseudo code for objective function claculation in documentation

### Gradient calculation of objective function

Here, first the gradient of the total cost function is calculated by expansion of \ref DerivOfObjectiveFunctionWRTModelParams "the derivative of the cost function with respect to the model parameters" and \ref CostFuncSummedOverAllTimePoints "the cost function summed over all time points"
The total gradient is given by

\f{align}{
\frac{\partial \mathcal{C}_\text{total}}{\partial \mathbf{R}} 
& = \sum_{t=1}^{N_i}
\frac{\partial \mathcal{C}_t}{\partial \mathbf{R}} \\
& = \sum_{t=1}^{N_i} \frac{\partial \mathbf{M}_t}{\partial\mathbf{R}}
\frac{\partial\mathcal{C}_t}{\partial\mathbf{M}_t}.
	\f}
Substituting \f$\frac{\partial\mathcal{C}_t}{\partial \mathbf{M}_t}\f$ with \ref ObjectiveFunction "the objective function" yields
\f{align}{
	\frac{\partial \mathcal{C}_\text{total}}{\partial \mathbf{R}} 
& = \sum_{t=1}^{N_i} 
\frac{\partial \mathbf{M}_t}{\partial\mathbf{R}}
\left( (1-\lambda)
       \frac{\partial \mathcal{S}_t}{\partial\mathbf{M}_t} + 
       \lambda
       \frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \right).
	\f}


Again, using the chain-rule to substitute \f$\frac{\partial \mathcal{S}_t}{\partial\mathbf{M}_t}\f$,
\f{align}{
	\frac{\partial \mathcal{C}_\text{total}}{\partial \mathbf{R}} 
& = \sum_{t=1}^{N_i} 
\frac{ \partial \mathbf{M}_t}{\partial\mathbf{R} }
\left( (1-\lambda) \frac{\partial \mathbf{I}_{T_t}}{\partial\mathbf{M}_t}\frac{\partial \mathcal{S}_t}{\partial\mathbf{I}_{T_t}} + \lambda
       \frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \right)\\
& = \sum_{t=1}^{N_i} 
\frac{ \partial \mathbf{M}_t}{\partial\mathbf{R} }
\left( 
(1-\lambda)
\frac{\partial \mathbf{I}_{T_t}}{\partial\mathbf{M}_t}
\frac{\partial \mathbf{P}_{\negthickspace A_t}}{\partial \mathbf{I}_{T_t}}
\frac{\partial \mathcal{S}_t}{\partial\mathbf{P}_{\negthickspace A_t}}  +
\lambda 
\frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \right)\\
& = \sum_{t=1}^{N_i} 
\frac{ \partial \mathbf{M}_t}{\partial\mathbf{R} }
\left( (1-\lambda) \frac{\partial \mathbf{I}_{T_t}}{\partial\mathbf{M}_t}
       A^*_t \left( \frac{\partial \mathcal{S}_t}{\partial\mathbf{P}_{\negthickspace A_t}} \right) 
+\lambda 
\frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \right).
\f}

In the next section the different elements of the derivative required for gradient-based optimisation are described in more detail.


#### Derivative: Regularisation with respect to the motion parameters

This term \f$\frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t}\f$ is well known from image registration and any existing implementation can be utilised here. 
Popular regularisation methods include for instance bending energy, linear elastic potential, and determinant of the Jacobian. 
 

#### Derivative: Similarity with respect to the simulated motion image

\anchor DerivSimWRTPartialImageData
\f[
\frac{\partial \mathcal{S}_t}{\partial\mathbf{P}_{\negthickspace A_t}}
=
\left[ 
\frac{\partial \mathcal{S}_t}{\partial p_{A_t,1}},
\ldots,
\frac{\partial \mathcal{S}_t}{\partial p_{A_t,N_p}}
\right]^\mathrm{T}
\f]
is a column vector of size \f$N_p \times 1\f$, where \f$N_p\f$ is the number of voxels in the simulated motion image \f$\mathbf{P}_{\negthickspace A_t}\f$. 
Each element of this vector is the gradient of the similarity measure with respect to that voxel (i.e. image intensity). 



If \f$\mathcal{S}\f$ measures similarity between the motion image and the simulated motion image in terms of the sum of squared differences as 
\f[
\mathcal{S}_\text{SSD} = \sum_{\mathbf{x\in\Omega}}\left(\mathbf{P}_{\negthickspace t}(\mathbf{x}) - \mathbf{P}_{\negthickspace A_t}(\mathbf{x})\right)^2,
\f]
then the derivative with respect to the simulated motion image intensity for each voxel is given by
\f[
 \frac{\partial \mathcal{S}}{\partial p_{A_t,y}}
    = -2\left(p_{A_t,y} - p_{t,y}\right).
\f]


#### Image acquisition and its adjoint 


\ref DerivSimWRTPartialImageData "The derivative of the image similarity with respect to the partial image data"  measures the gradient direction in the space of the raw 
(or partial, or unsorted) image data. In order transform this into the space of the full, reconstructed image space, the adjoint of the image acquisition
function \f$A\f$ may be used. 

The image acquisition process can be described by a time-dependent function \f$A_t\f$, which maps the unknown underlying image \f$\mathbf{I}_t\f$ -- or the 
transformed reference state image \f$\mathbf{I}_{T_t}\f$ into the acquisition space. 
\anchor simulationOfImageAcquisition
\f[
\mathbf{P}_{\negthickspace t} = A_t(\mathbf{I}_t) + \varepsilon_t
\f]
As a consequence, the image acquisition \f$A_t\f$ can be expressed as an \f$N_p \times N_v\f$ matrix
\f[
A_t = 
\begin{pmatrix}
a_{t,1,1} & \ldots & a_{t,1,N_v} \\
\vdots & \ddots & \vdots \\
a_{t,N_p,1} & \ldots & a_{t,N_p,N_v}
\end{pmatrix}
\f]
where each element \f$a_{t,i,j}\f$ is the contribution of voxel \f$j\f$ in the transformed reference state image to the partial image data voxel at position \f$i\f$.
\f[
a_{t,i,j} = \frac{\partial p_{A_t,i}}{\partial i_{T_t,j}}
\f]

Using the hermitian adjoint, i.e. complex conjugate, \f$A^*_t = \left(\bar{A_t}\right)^\mathrm{T}\f$, the similarity difference measured in the 
partial image space can be mapped back to full image space. Consequently, \f$A_t^*\f$ is of size \f$N_v\times N_p\f$ when written in matrix form.

\anchor adjointAsDerivative
\f[
A_t^* = \frac{\partial \mathbf{P}_{\negthickspace A_t}}{\partial \mathbf{I}_{T_t}}
\f]
Note, appropriate masking is required in the implementation of this mapping, depending on the reach of the partial image data on the full image data. 





#### Derivative: Transformed reference state image with respect to motion parameters

This derivative , \f$\frac{\partial \mathbf{I}_{T_t}}{\partial\mathbf{M}_t}\f$ gives the changes of the intensity of the transformed 
reference image with respect to the motion parameters and is thus dependent on the transformation model used. In the special case of 
non-parametric registration, this is given by the warped spatial image gradient. 

Hence, the derivative given in matrix notation is of size \f$N_m \times N_v\f$.
\f[
\frac{\partial \mathbf{I}_{T_t}}{\partial \mathbf{M}_t} 
= \begin{pmatrix}
\frac{\partial i_{T_t,1}}{\partial m_{t,1}}
& \ldots & \frac{\partial i_{T_t,N_v}}{\partial m_{t,1}}\\
\vdots & \ddots & \vdots \\
\frac{\partial i_{T_t,1}}{\partial m_{t,N_m}} & \ldots & 
\frac{\partial i_{T_t,N_v}}{\partial m_{t,N_m}}
\end{pmatrix}
\f]
Note that \f$N_m\f$ is the number of motion parameters which, in case of a parametric transformation is the product of the number of 
motion parameters (e.g. control points for a B-spline transformation), and in case of a non-parametric transformation the number of 
deformation vectors times the number of dimensions. 

In case of a parametric transformation, such as a B-spline transformation, the chain rule is applied once more
\anchor imageGradWRTMotionParameters
\f[
\frac{\partial \mathbf{I_{T_t}}}{\partial \mathbf{M}_t} 
=
\frac{\partial \textbf{DVF}_t}{\partial \mathbf{M_t}}
\frac{\partial \mathbf{I}_{T_t}}{\partial \textbf{DVF}_t}.
\f]
which means for each matrix entry
\f[
\frac{\partial i_{T_t,j}}{\partial m_{t,i}}
= 
\frac{\partial i_{T_t,j}}{\partial \textit{dvf}_{t,x}}
\frac{\partial \textit{dvf}_{t,x}}{\partial m_{t,i}}.
\f]
This is the multiplication of the warped image gradient with the contribution of each control point to the transformation at a given position 
\f$x\f$. The latter one can be calculated by convolving the deformation vector field with the B-spline base function and sampling accordingly. 


#### Derivative: Motion parameters wirth respect ot model parameters

The \ref correspondenceModelQuadratic "quadratic correspondence model" can be generalised for any polynomial function as 
\anchor correspondenceModelPolynomial
\f[
\mathbf{M}_t = \sum_{i=1}^{N_s} \mathbf{R}_i s_{i,t}
\f]
with the corresponding partial derivative with respect to the model parameters
\anchor partialMotionWRTModel
\f[
\frac{\partial\mathbf{M}_t}{\partial\mathbf{R}_i} = s_{i,t}.
\f]
In matrix notation \ref partialMotionWRTModel "the above equation" can be written as a vertical stack of diagonal matrices:
\f[
\frac{\partial \mathbf{M}_t}{\partial\mathbf{R}}
=\begin{pmatrix}
\frac{\partial\mathbf{M}_t}{\partial\mathbf{R}_1}\\
\frac{\partial\mathbf{M}_t}{\partial\mathbf{R}_2}\\
\vdots\\
\frac{\partial\mathbf{M}_t}{\partial\mathbf{R}_{N_s}}
\end{pmatrix}
=
\begin{pmatrix}
s_{1,t} \\
&\ddots\\
&&s_{1,t}\\
&\vdots\\
s_{N_{s},t} \\
&\ddots\\
&&s_{N_{s},t}\\
\end{pmatrix}
\f]

This matrix is of size \f$N_r \times N_m\f$ and completes the collection of partial derivatives required for the gradient calculation. 
A summary of the complete procedure how to compute the gradient is presented in listing ...

\todo Include pseudo code for objective function gradient calculation. 











\tableofcontents
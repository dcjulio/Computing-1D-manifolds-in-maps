# Computing-1D-manifolds-in-maps

LICENSE: In case you use GrowFundCurv1D in your own research, please, give credit by citing [1]. 

The package GrowFundCurv1D is an implementation in Matlab of an algorithm for the numerical computation of one-dimensional (un)stable manifolds and intersection points in maps. This is joint work with Hinke. M. Osinga and Bernd Krauskopf. 

The package GrowFundCurv1D is an implementation in Matlab of the algorithm presented in [1] for the numerical computation of one-dimensional (un)stable manifolds and intersection points in maps, based on the manifold growth algorithm from [2]. This manual serves as a practical guide for utilising the package GrowFundCurv1D without delving into extensive algorithmic details. We refer to [1] for a comprehensive discussion of the algorithm and its accuracy constraints, as well as more details on its performance.

The package GrowFundCurv1D has been designed for three-dimensional maps, and uses an initial short segment of a stable or unstable manifold associated with a fixed or periodic point. This initial segment can be obtained either from the linear approximation, given by the eigenvector associated with the respective stable or unstable eigenvalue, or from a previous continuation. The algorithm first identifies the last fundamental domain on the initial segment, based on which the iterative process of growing the manifold starts. As a post-processing step, the algorithm computes the intersection points of the computed manifold with a pre-specified plane as an ordered set.

The package GrowFundCurv1D comprises a series of routines, explained in the comprehensive step-by-step manual available as GrowFundCurv1D_manual.pdf. The package and routines are available in the folder GrowFundCurv1D_code/. The folder contains the Matlab script file GrowFundCurv1D_demo.m that demonstrates the algorithm with a specific example, the data used to obtain the first fundamental domain of the manifold, and the folder GrowFundCurv1D_functions/ with required functions.

The demo GrowFundCurv1D_demo.m, was tested using Matlab [version 9.12 (R2022a)]. This example computes a one-dimensional stable manifold of a fixed point for a three-dimensional Hénon-like map, as defined in [1]. Note that, with appropriate changes to the accuracy settings, the algorithm can accurately compute manifolds not only for the fixed points of the map itself, but also for up to its fourth iterate, without losing resolution.

[1] D. C’Julio, B. Krauskopf, and H.M. Osinga. Computing parametrised large intersection sets of 1D invariant manifolds: a tool for blender detection. Preprint available from https://www.math.auckland.ac.nz/~hinke/preprints/cko_algorithm.html, 2023.

[2] B. Krauskopf and H. M. Osinga. Growing 1D and quasi-2D unstable manifolds of maps. Journal of Computational Physics, 146(1): 404–419, 1998.

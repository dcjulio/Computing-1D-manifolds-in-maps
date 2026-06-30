# Computing-1D-manifolds-in-maps


NOTE:<br> This is version v3.0.1.<br>
See the file CHANGELOG for changes between sub-versions.<br>
To browse other previous versions, click on "Releases" on the right side of the page. There, all the published releases are shown in reverse chronological order.

---
LICENSE: In case you use GrowFundCurv1D in your own research, please, give credit by citing [1]. If you use the capability of the algorithm to extract a pseudo orbit, please cite [2].

The package GrowFundCurv1D is a Matlab implementation of the computational tool described in [1], based on the manifold growth algorithm from [3] for computing one-dimensional (un)stable manifolds of fixed points in three-dimensional maps. It identifies the intersection points of the manifold with a section, stores the indices of approximate preimages for all mesh points, enabling the reconstruction of pseudo-orbits, and visualises the results. The package GrowFundCurv1D is joint work with Hinke. M. Osinga and Bernd Krauskopf. Version 3.0 is also joint work with Sanaz Amani and Sam Doak.

The algorithm begins by computing an initial segment of the manifold via linear interpolation in the direction of the eigenvector associated with the fixed point, using a small step delta. Starting from this segment, the algorithm identifies a fundamental domain and initiates the iterative growth process. During this process, the algorithm stores the approximate preimages of each mesh point to enable a pseudo-orbit construction.
As a post-processing step, the intersection points of the computed manifold with a pre-specified plane are determined as an ordered set. For a comprehensive explanation of the algorithm's performance and accuracy constraints, see [1]. For documentation on pseudo-orbit computation and its application to boundary value problems (BVP) and parameter continuation, please refer to [2].

The package GrowFundCurv1D comprises a series of routines, explained in the manual available as GrowFundCurv1D_manual.pdf. The package and routines are available in the folder GrowFundCurv1D_code/. The folder contains the Matlab script file GrowFundCurv1D_demo.m that demonstrates the algorithm with a specific example, and the folder GrowFundCurv1D_functions/ with the required functions.

The demo GrowFundCurv1D_demo.m, was tested using Matlab [from version 9.12 (R2022a) to version 24.2 (R2025b)]. This example computes a one-dimensional stable manifold of a fixed point for a three-dimensional Hénon-like map, as defined in [1]. Note that, with appropriate changes to the accuracy settings, the algorithm can accurately compute manifolds not only for the fixed points of the map itself, but also for up to its fourth iterate, without losing resolution.

[1] D. C’Julio, B. Krauskopf, and H.M. Osinga. Computing parametrised large intersection sets of 1D invariant manifolds: a tool for blender detection. Numerical Algorithms 96(3): 1079–1108, 2024.

[2] D. C'Julio, S. Amani, S. Doak, B. Krauskopf, H.M. Osinga. Constructing pseudo-orbits on invariant manifolds of maps as seeds for boundary value problems
Preprint, University of Auckland (2025)

[3] B. Krauskopf and H. M. Osinga. Growing 1D and quasi-2D unstable manifolds of maps. Journal of Computational Physics, 146(1): 404–419, 1998.

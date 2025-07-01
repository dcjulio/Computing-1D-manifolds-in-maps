### What’s new in v2.0.0
1. Manifold initialisation is now handled directly in the GrowFundCurv1D function — it no longer requires an initial segment computed externally.
2. Updated syntax for accessing manifold points. They are now manif.points.pos and manif.points.neg, replacing the old manif.pointspos and manif.pointsneg.
3. Users can define both the final arclength of the manifold to compute and/or a maximum number of fundamental domain iterations. The algorithm will grow the manifold up to the specified arclength unless it reaches the iteration limit first.
4. You can now add new branches to an already computed manifold by calling the function "add_branch".
5. For orientation-reversing manifolds (with negative eigenvalues), the algorithm now computes both branches using the first iterate of the map. In v1.0.0, it used the second iteration of the map and handled each branch separately.
6. The algorithm now stores the indices corresponding to each fundamental domain or segment.
7. The inter_planes function no longer adds intersection points directly to the manifold. Instead, it records each intersection point along with the index of the previous point in the manifold that intersected the plane.

### New Warnings:
<b>Discontinuity Alert:</b> A WARNING will now be printed if the distance between the last point of the previous segment and the first point of the new segment exceeds Delta_min. This indicates a potential discontinuity in the manifold.

<b>Accuracy Check:</b> The system will now issue a WARNING if the accuracy conditions are not met between two concatenated fundamental domains, signalling a potential precision issue.

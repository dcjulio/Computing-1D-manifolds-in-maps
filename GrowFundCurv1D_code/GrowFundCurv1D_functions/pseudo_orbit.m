function orbit=pseudo_orbit(manif, idx, branch)
% idx: the index to compute the pseudo orbit with
% branch: it can be 'pos' or 'neg'

    i=0;

    while idx~=0
        i=i+1;

        %where the point comes from
        orbit.name{i} = [manif.fixp.name '_' branch];
        %coordinates of the point
        orbit.x(i) = manif.points.(branch).x(idx);
        orbit.y(i) = manif.points.(branch).y(idx);
        orbit.z(i) = manif.points.(branch).z(idx);
        %index of the point in the manifold
        orbit.idx(i) = idx;

        % information for the next preimage
        idx = manif.points.(branch).pre_idx(idx); %index of the manifold preimage
        branch   = manif.points.(branch).pre_branch;         %branch of the preimage
        
    end

end
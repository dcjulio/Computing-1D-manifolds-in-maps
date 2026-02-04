function orbit=eps_pseudo_orbit(manif, idxpoint, branch)
% idxpoint: the index to compute the pseudo orbit with
% branch: it can be 'pos' or 'neg'

    i=0;

    while idxpoint~=0
        i=i+1;

        %where the point comes from
        orbit.name{i} = [manif.fixp.name '_' branch];
        %coordinates of the point
        orbit.x(i) = manif.points.(branch).x(idxpoint);
        orbit.y(i) = manif.points.(branch).y(idxpoint);
        orbit.z(i) = manif.points.(branch).z(idxpoint);
        %index of the point in the manifold
        orbit.idxpoint(i) = idxpoint;

        % information for the next preimage
        idxpoint = manif.points.(branch).pre_idx(idxpoint); %index of the manifold preimage
        branch   = manif.points.(branch).pre_branch;         %branch of the preimage
        
    end

end
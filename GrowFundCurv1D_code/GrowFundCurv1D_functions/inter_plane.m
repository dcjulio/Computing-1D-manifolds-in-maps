function manif=inter_plane(manif,angle)
%               (as input)
%               manif.inter.angle    => the angle of the plane %0 is y==0. Between 0 and +-pi. (pi/2: x==0 (y>0))
%               manif.points         => points.x points.y points.z

%               (as output)
%               .vec      => the vector that generates the plane
%               .plane    => the equation of the plane
%               .idx and/or .idxneg, .idxpos => idx if is orientable, 
%                                               idxneg(pos) or non-orientable

format long

npoint=3; %number of points to consider to find intersect with spline
manif.inter.angle=angle;
angle = manif.inter.angle; %angle of the plane %0 is y==0. Between 0 and +-pi. (pi/2: x==0 (y>0))
vect=[cos(angle),sin(angle)]; %vector for the plane
    
[n, d] = numden(sym(manif.inter.angle/pi));
if n==0
    manif.inter.angle=sprintf('%ipi',n);
else
    manif.inter.angle=sprintf('%i/%ipi',n,d);
end

%% just to save the name of the plane (for example '{x=2y}')
if angle>0
    side=', y>0';
elseif angle<0
	side=', y<0';
end

if mod(angle,pi)==pi/2 %{x=0}
    manif.inter.plane=sprintf('{x=0%s}',side); 
elseif angle==0
	manif.inter.plane='{y=0}';
else
    coefficients = polyfit([0, vect(1)], [0, vect(2)], 1); 
    a = coefficients(1);
    
    if round(a*1e2)==100
        manif.inter.plane=sprintf('{y=x%s}',side); 
    elseif round(a*1e2)==-100
        manif.inter.plane=sprintf('{y=-x%s}',side); 
    else
        manif.inter.plane=sprintf('{y=%.2fx%s}',a,side); 
    end
end


%% intersection with plane

if isfield(manif.points,'pos')
    [manif.inter.points.pos]=plane_inter(manif.points.pos, angle, npoint);
end
if isfield(manif.points,'neg')
    [manif.inter.points.neg]=plane_inter(manif.points.neg, angle, npoint);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%
function [inter_points]=plane_inter(points, angle, npoint)

    vec=[cos(angle),sin(angle)]; %vector for the plane
    p0=[0,0]; %point in the plane
    N=[vec(2),-vec(1)]; %normal vector

    % points n and n+1 (to later check when they cross the plane)
    % we only check x and y because the plane doesnt depend on z
    X1=points.x(1:end-1);
    Y1=points.y(1:end-1);
    X2=points.x(2:end);
    Y2=points.y(2:end);
    
    % dot((x,y)-p0,n). %if (x,y) is in the plane then eq==0
    plane_eq =  @(x,y) (x-p0(1)).*N(1) + (y-p0(2)).*N(2); 
%% Finding when the points cross the plane
%----------- Finding the indices of points that are crossing the plane, in any direction
    % n > plane   and    n+1 < plane
    idx1 = intersect(find(plane_eq(X1,Y1) >= 0), find(plane_eq(X2,Y2) < 0));
    % n < plane   and    n+1 > plane
    idx2 = intersect(find(plane_eq(X1,Y1) <= 0), find(plane_eq(X2,Y2) > 0));
    
    idx_all=union(idx1,idx2); %joining those indices
    idx_all=idx_all(idx_all<numel(X1));
    
%---%------- Find only the index of points that are crossing the correct side of the plane (because is a half plane, when y>0 or y<0)

%---%------- if the plane is not y==0
    if vec(2)~=0 
%---%------- 
        if angle>0
            idx.all=intersect(idx_all,find(Y1>=0));
        elseif angle<0
            idx.all=intersect(idx_all,find(Y1<=0));
        end
%---%------- if the plane is y==0
    else
%---%------- 
        idx.all=sort(idx_all);
%---%------- 
    end
%---%-------

% they have to have at least "npoint" after the index
idx.all = idx.all(idx.all+npoint<numel(points.x));
idx.all = idx.all(idx.all>npoint); %and "npoint" before the index

%% Finding when the points cross the plane

% we will save the coordinates of the intersection points herer
inter_points.x=zeros(1,numel(idx.all));
inter_points.y=zeros(1,numel(idx.all));
inter_points.z=zeros(1,numel(idx.all));

inter_points.idx=idx.all(idx.all>1); %index before crossing the plane

%---%------- If the manif cross the plane at least one time, we compute the intersection coordinates
    if numel(idx.all)>=1 
%---%-------     

        %define the segment of manifold that is passing through the plane
        X=zeros(numel(idx.all),2*npoint); 
        Y=zeros(numel(idx.all),2*npoint);
        Z=zeros(numel(idx.all),2*npoint);

%---%---%-- npoints on one side of the plance, npoints on the other side
        for k=1:(2*npoint)
            X(:,k)=points.x(idx.all+npoint-(k-1))';
            Y(:,k)=points.y(idx.all+npoint-(k-1))';
            Z(:,k)=points.z(idx.all+npoint-(k-1))';
        end
%---%---%-- 

        % defining the structure with points on the plane
        inter= struct('x',num2cell(X,2),'y',num2cell(Y,2),'z',num2cell(Z,2));

%---%---%-- define the spline (makima) for each segment
        t=1:2*npoint;
        for kk=1:numel(inter) 
%---%---%--

            %using makima to interpolate the points in the plane for each coordinate
            pp_x = makima(t,inter(kk).x);
            pp_y = makima(t,inter(kk).y);
            pp_z = makima(t,inter(kk).z); 

            %coeficients
            c_x=pp_x.coefs(npoint,:);
            c_y=pp_y.coefs(npoint,:);
            c_z=pp_z.coefs(npoint,:);

            %breaks
            breaks_x = pp_x.breaks(npoint);
            breaks_y = pp_y.breaks(npoint);
            breaks_z = pp_z.breaks(npoint);

%---%---%--% Defining the equations to find the points on the plane
%---%---%--% If we have more than one point per side of the plane then is not linear
            if npoint>1
%---%---%--% 
                % ecuation of spline
                XX = @(T) c_x(:,1).*(T-breaks_x').^3+c_x(:,2).*(T-breaks_x').^2+c_x(:,3).*(T-breaks_x')+c_x(:,4);
                YY = @(T) c_y(:,1).*(T-breaks_y').^3+c_y(:,2).*(T-breaks_y').^2+c_y(:,3).*(T-breaks_y')+c_y(:,4);
                ZZ = @(T) c_z(:,1).*(T-breaks_z').^3+c_z(:,2).*(T-breaks_z').^2+c_z(:,3).*(T-breaks_z')+c_z(:,4);
%---%---%--% with only two point (one each side of plane) we can only use a linear equation
            else %npoint=1
%---%---%--% 
                % ecuation linear
                XX = @(T) c_x(:,1).*(T-breaks_x')+c_x(:,2);
                YY = @(T) c_y(:,1).*(T-breaks_y')+c_y(:,2);
                ZZ = @(T) c_z(:,1).*(T-breaks_z')+c_z(:,2);
%---%---%--% 
            end
%---%---%--% 

            % Using fsolve to find the points on the plane
            options = optimoptions(@fsolve,'FunctionTolerance',1e-18,'OptimalityTolerance',1e-18,'StepTolerance',1e-18,'Display','off');
            tt = fsolve(@(T) plane_eq(XX(T),YY(T)) , npoint, options);

            %this are the new points on the plane
            inter_points.x(kk)=XX(tt);
            inter_points.y(kk)=YY(tt);
            inter_points.z(kk)=ZZ(tt);

%---%---%-- 
        end
%---%---%-- 
    end

end %function plane_inter

end %function plane_angle_all_manif_makima
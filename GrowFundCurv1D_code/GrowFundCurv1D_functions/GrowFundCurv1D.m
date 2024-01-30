function Manif=GrowFundCurv1D(manif,opts)
    
    tic

    %%
    %saves the information of the system stored in the structure 'manif' and 'opts'
    inf_sys=manif.inf_sys;
    thesystem=opts.thesystem;

    %% Warnings in case the specifications are not correct for the computation of this manifold
    if strcmp(manif.orientability,'orientation-reversing')
        if mod(manif.growinf.mapiter,2)==1
            fprintf('\n----  Warning!  Manifold is non-orientable and the number of iterations of the map has to be even\n');
            manif.growinf.mapiter=manif.growinf.mapiter+1;
            fprintf('----  Number of total iterations of the map has been updated to %i ---- \n\n',manif.growinf.mapiter);
            prompt = "... Press any key to continue... \n\n";
            x = input(prompt);
        end
    end

    if floor(manif.growinf.funditer/manif.growinf.mapiter)~=manif.growinf.funditer/manif.growinf.mapiter
        fprintf('\n----  Warning! Number of total iteration of the manifold has to be multiply of iterations of the map (%i) ---- \n',manif.growinf.mapiter);
        manif.growinf.funditer=max(manif.growinf.mapiter,floor(manif.growinf.funditer/manif.growinf.mapiter)*manif.growinf.mapiter);
        fprintf('----  Number of total iterations of the manifold has been updated to %i ---- \n\n',manif.growinf.funditer);
        prompt = "... Press any key to continue... \n\n";
        x = input(prompt);
    end

%% Printing the general information of this run

    names=fieldnames(inf_sys.par);
    
    fprintf('\n----');
    for k=1:length(names)
        fprintf(' %s:%0.2f ',names{k},inf_sys.par.(names{k}));
    end
    fprintf('----');

    fprintf('\n----  Total iterations of the fundamental domain: %i ---- ',manif.growinf.funditer);
    fprintf('\n----  Iterations of map used: %i ---- ',manif.growinf.mapiter);
    fprintf('\n\n')
    fprintf('----- Acc. Conditions ----- ')
    fprintf('\n| AlphaMax:%.2e',manif.growinf.alphamax)
    fprintf('\n| DeltaAlphaMax:%.2e ',manif.growinf.deltalphamax)
    fprintf('\n| DeltaMin:%.2e ' ,manif.growinf.deltamin)
    fprintf('\n| DeltaMax:%.2e \n' ,manif.growinf.deltamax)
    fprintf('---------- \n')
    
    
    %-------------------------------------------------------------------
    %-------------------------------------------------------------------
    % initializing the fields for the initial information
    manif.growinf.runinf.rem_deltamin=0; %how many points are removed because of deltamin
    manif.growinf.runinf.rem_nan=0; %how many points are removed because of NaN values
    manif.growinf.runinf.rem_inf=0; %how many points are removed at infinity because of duplication
    manif.growinf.runinf.add_alphamax=0; %how many points are added because of alpha_max
    manif.growinf.runinf.add_deltamax=0;%how many points are added because of delta_max
    manif.growinf.runinf.add_deltalphamax=0;  %how many points are added because of (delta alpha)_max
    manif.growinf.runinf.npoints_initial_final=numel(manif.points.x); %how many points in the manifold at the start and at the end
    manif.growinf.runinf.arc_initial_final=manif.points.arc(end); %arclength at the start and at the end



    %find the last fundamental domain for matcont manifold
    fund_end.x=manif.points.x(end);
    fund_end.y=manif.points.y(end);
    fund_end.z=manif.points.z(end);
    if strcmp(manif.stab,'Smanifold')
        init_fund = thesystem.mapping(fund_end,'Umanifold',opts);
    elseif strcmp(manif.stab,'Umanifold')
        init_fund = thesystem.mapping(fund_end,'Smanifold',opts);
    end
    fund_domain.x=[init_fund.x fund_end.x];
    fund_domain.y=[init_fund.y fund_end.y];
    fund_domain.z=[init_fund.z fund_end.z];

% find the closesst point to the pre-image of the last point
dist=(manif.points.x-init_fund.x).^2+(manif.points.y-init_fund.y).^2+(manif.points.z-init_fund.z).^2;
[~,ind]=min(dist);

%where to add the extra point (I only check x coordinate because how my map behaves
if (init_fund.x-manif.points.x(ind)>0 && init_fund.x-manif.points.x(ind+1)<0) || (init_fund.x-manif.points.x(ind)<0 && init_fund.x-manif.points.x(ind+1)>0)
    manif.points.x=[manif.points.x(1:ind) init_fund.x manif.points.x(ind+1:end)];
    manif.points.y=[manif.points.y(1:ind) init_fund.y manif.points.y(ind+1:end)];
    manif.points.z=[manif.points.z(1:ind) init_fund.z manif.points.z(ind+1:end)];
else
   	manif.points.x=[manif.points.x(1:ind-1) init_fund.x manif.points.x(ind:end)];
    manif.points.y=[manif.points.y(1:ind-1) init_fund.y manif.points.y(ind:end)];
    manif.points.z=[manif.points.z(1:ind-1) init_fund.z manif.points.z(ind:end)];   
end

%%

Manif=manif;
Manif.points.arc = arclength(Manif.points);
fprintf('\n Starting arclength  %f\n',Manif.points.arc(end))

% Here we store the computation of the fundamental domains
idx_fund(1)=intersect(find(Manif.points.x==fund_domain.x(end-1)),find(Manif.points.y==fund_domain.y(end-1)));
idx_fund(2)=intersect(find(Manif.points.x==fund_domain.x(end)),find(Manif.points.y==fund_domain.y(end)));


%define the last fundamental domain
fund.points.x=Manif.points.x(idx_fund(1):idx_fund(2));
fund.points.y=Manif.points.y(idx_fund(1):idx_fund(2));
fund.points.z=Manif.points.z(idx_fund(1):idx_fund(2));
fund.points.arc = arclength(fund.points);

fprintf(' Starting fundamental domain arclength  %f\n',fund.points.arc(end))

%----------- Starting the loop
for iter=1:floor(manif.growinf.funditer/manif.growinf.mapiter) %how many times do we have to iterate the algorithm
    
    %mapping the points
    mappoints = thesystem.mapping(fund.points,manif.stab,opts);
    
    %% STARTING THE ALGORITHM

    %----------- Removing points at infinity or NaN values

    %Delete NaN values
    nan_idx=union(union(find(isnan(mappoints.x)),find(isnan(mappoints.y))),find(isnan(mappoints.z)));
    Manif.growinf.runinf.rem_nan=Manif.growinf.runinf.rem_nan+numel(nan_idx);
    
    mappoints.x(nan_idx)=[]; 
    mappoints.y(nan_idx)=[]; 
    mappoints.z(nan_idx)=[]; 
    
    fund.points.x(nan_idx)=[]; 
    fund.points.y(nan_idx)=[]; 
    fund.points.z(nan_idx)=[]; 
    
    
    %Delete points at infinity
    inf_idx=union(find(sqrt(mappoints.x.^2+mappoints.y.^2)==1),find(abs(mappoints.z)==1));
    Manif.growinf.runinf.rem_inf=Manif.growinf.runinf.rem_inf+numel(inf_idx);
    
    mappoints.x(inf_idx)=[]; 
    mappoints.y(inf_idx)=[]; 
    mappoints.z(inf_idx)=[]; 
    
    fund.points.x(inf_idx)=[]; 
    fund.points.y(inf_idx)=[]; 
    fund.points.z(inf_idx)=[]; 


    fprintf('\n ITERATION OF THE FUNDAMENTAL DOMAIN %i\n',iter*manif.growinf.mapiter)

    
%---%----------- Adding points depending on Acc. Cond.  
    
    %initializing the structures to add points
    add_acc=struct(); 
    newpoints=struct();
    mapnewpoints=struct();
    
    
	add_acc.iter=0;
    add_acc.failed=[]; % points that failed acc cond in last loop
    add_acc.loop=true; % still doing the while loop 


    % Interpolation mesh
    fund_initial=fund.points; %starting mesh of the fundamental domain
    t_initial=0:1/(numel(fund_initial.x)-1):1; % parametrization for meshpoints

%---%-------------- Loop of the same mesh checking acc cond (this adds points)
    while add_acc.loop 
%---%--------------
        add_acc.loop=false; %to stop while loop % if at least one point is added this turns true
        add_acc.iter=add_acc.iter+1;

        % Interpolating points from previous fundamental domain.
         if add_acc.iter==1
            tt=t_initial;
            t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of interpolated points
            interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
            mapinterp = thesystem.mapping(interp,manif.stab,opts); % interpolated image
         else
             tt=sort([tt t_interp(add_acc.add)]); %parametrization of (new) mesh points
             t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of (new) interpolated points
             interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
             mapinterp = thesystem.mapping(interp,manif.stab,opts); % interpolated image
         end
%%
        add_acc.add=[]; %points we are going to add

        % idx of points to check acc cond
        if add_acc.iter==1
            for_idx=2:(numel(mappoints.x)-1);
        else
            for_idx=add_acc.failed;
        end      

        fprintf('  loop number %i (points to check %i...',add_acc.iter,numel(for_idx));

        
%-------%---------- Going through the points that failed
        million=0;
        for k=for_idx
%-------%----------      
            % a flag for when # million points have been checked
            million=million+1;

            if floor(million/1000000)==ceil(million/1000000) %is integer?
                fprintf(' -checkpoint %i million points checked-...',floor(million/1000000));
            end

                
            % coordinates of mapped points
            add_acc.p0=[mappoints.x(k-1), mappoints.y(k-1), mappoints.z(k-1)];
            add_acc.p1=[mappoints.x(k), mappoints.y(k), mappoints.z(k)]; % the point we are actually looking at
            add_acc.p2=[mappoints.x(k+1), mappoints.y(k+1), mappoints.z(k+1)];
            
            % Distance btw points
            add_acc.delta0=norm(add_acc.p1-add_acc.p0); % before
            add_acc.delta2=norm(add_acc.p1-add_acc.p2); % after
            add_acc.alpha = angles(add_acc.p0,add_acc.p1,add_acc.p2); % angle btw points
            
            % points btw p0p1 and p1p2 in the interpolated points
            add_acc.p0_new=[mapinterp.x(k-1), mapinterp.y(k-1), mapinterp.z(k-1)];
            add_acc.p2_new=[mapinterp.x(k), mapinterp.y(k), mapinterp.z(k)];
    

%-----------%------ Adding points


%-----------%------ If delta > deltamax 

            %------ If it is the second point in the mesh
            %------ Check the first delta and add a point before
            if k==2 && add_acc.delta0>manif.growinf.deltamax  

                %add point p01
                add_acc.add    =[add_acc.add k-1]; %idx of the point we are going to add
                add_acc.loop=true; % We have to check if we need to put more points in the mesh
                Manif.growinf.runinf.add_deltamax=Manif.growinf.runinf.add_deltamax+1;
                
            %------ Check the second delta and add a point after 
            elseif add_acc.delta2>manif.growinf.deltamax %k>2 && add_acc.delta2>manif.growinf.deltamax
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.growinf.runinf.add_deltamax=Manif.growinf.runinf.add_deltamax+1;
           
                
                
%-----------%------ If alpha > alphamax  or   Delta*alpha > Delta*alpha max

            %------ If alpha fails or BOTH Delta*alpha fail
            %------ Choose where to add a point.
            elseif add_acc.alpha>=manif.growinf.alphamax || (add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax) %1 

            %------ Add point where Delta>Deltamin
                
                %-- If only Delta0>Deltamin
                %-- Add a point btw p0 and p1
                if add_acc.delta0>manif.growinf.deltamin && add_acc.delta2<manif.growinf.deltamin 
                    % Add a point if we didnt added the point in the previous acc cond checks
                    if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                        add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax) 
                            Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                        else
                            Manif.growinf.runinf.add_alphamax=Manif.growinf.runinf.add_alphamax+1;
                        end
                    end
                end
                
                %-- If only Delta2>Deltamin
                %-- Add a point btw p1 and p2
                if add_acc.delta2>manif.growinf.deltamin && add_acc.delta0<manif.growinf.deltamin
                    %add point p12
                    add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    if (add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax) 
                        Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                    else
                        Manif.growinf.runinf.add_alphamax=Manif.growinf.runinf.add_alphamax+1;
                    end
                end
                
                %-- If both Delta0 and Delta 2 > Deltamin
                %-- Choose where to add point
                if add_acc.delta2>manif.growinf.deltamin && add_acc.delta0>manif.growinf.deltamin
                    add_acc.alpha0_new = angles(add_acc.p0_new,add_acc.p1,add_acc.p2); % angle btw points
                    add_acc.alpha2_new = angles(add_acc.p0,add_acc.p1,add_acc.p2_new); 

                    if add_acc.alpha0_new < add_acc.alpha2_new 
                        % Add a point if we didnt added the point in the previous acc cond checks
                        if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                            add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                            add_acc.loop=true; %idx of the point we are going to add
                            if (add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax) 
                                Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                            else
                                Manif.growinf.runinf.add_alphamax=Manif.growinf.runinf.add_alphamax+1;
                            end
                        end
                    else
                        %add point p12
                        add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax) 
                            Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                        else
                            Manif.growinf.runinf.add_alphamax=Manif.growinf.runinf.add_alphamax+1;
                        end
                    end
                end
                
                
                
%-----------%------ If just one Delta*alpha > Delta*alpha max ( and alpha < alphamax (previous elseif is when alpha > alphamax )

            %------ Delta0*alpha > max, and Delta0 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta0*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta0>manif.growinf.deltamin
                % Add a point if we didnt added the point in the previous acc cond checks
                if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                    add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                end
                
            %------ Delta2*alpha > max, and Delta2 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta2*add_acc.alpha>=manif.growinf.deltalphamax && add_acc.delta2>manif.growinf.deltamin
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.growinf.runinf.add_deltalphamax=Manif.growinf.runinf.add_deltalphamax+1;
                
                
%-----------%------              
            end    % (if loop) Adding points     
%-----------%------      
%-------%---------- 
        end       % (for loop) Going through the points that failed 
%-------%----------

        newpoints.x=interp.x(add_acc.add);
        newpoints.y=interp.y(add_acc.add);
        newpoints.z=interp.z(add_acc.add);
        
        mapnewpoints.x=mapinterp.x(add_acc.add);
        mapnewpoints.y=mapinterp.y(add_acc.add);
        mapnewpoints.z=mapinterp.z(add_acc.add);
        

       fprintf(' added points: %i) \n',numel(add_acc.add));

        
        % Add the points 
        if ~isempty(add_acc.add)
            
            % get updated idx of failed points
            plus=0:(numel(add_acc.add)-1);
            add_acc.failed = unique(sort([add_acc.add+plus add_acc.add+plus+1 add_acc.add+plus+2]));
            add_acc.failed=add_acc.failed(add_acc.failed>1);
            add_acc.failed=add_acc.failed(add_acc.failed<numel(mappoints.x)+numel(add_acc.add));
            


            % add points in the mapped manifold and in the old manifold
            mappoints.x=insert(mappoints.x,mapnewpoints.x,add_acc.add);
            mappoints.y=insert(mappoints.y,mapnewpoints.y,add_acc.add);
            mappoints.z=insert(mappoints.z,mapnewpoints.z,add_acc.add); 

            fund.points.x=insert(fund.points.x,newpoints.x,add_acc.add);
            fund.points.y=insert(fund.points.y,newpoints.y,add_acc.add);
            fund.points.z=insert(fund.points.z,newpoints.z,add_acc.add); 

        end
        
%---%--------------       
    end           % (while loop) Checking acc cond 
%---%--------------

%---%----------- Final section: save info

    fund.points.x=mappoints.x;
    fund.points.y=mappoints.y;
    fund.points.z=mappoints.z;
    fund.points.arc = arclength(fund.points);
    
    %add new branch to the entire manifold
    Manif.points.x=[Manif.points.x fund.points.x(2:end)];
    Manif.points.y=[Manif.points.y fund.points.y(2:end)];
    Manif.points.z=[Manif.points.z fund.points.z(2:end)];
    Manif.points.arc = arclength(Manif.points);
    
    
%---%----------- END Final section: save info
end
%--------------- END adding points



%--------------- Removing points depending on Deltamin
% We remove points that are too close together. They are a result of the points mapped, not the algorithm that adds points
% We do this at the end because they maybe were useful for interpolating points

rem_deltamin=struct(); % structure for removing points

% find distance between points and save index of points < dmin
rem_deltamin.delta=vecnorm([Manif.points.x(1:end-1); Manif.points.y(1:end-1); Manif.points.z(1:end-1)] - [Manif.points.x(2:end); Manif.points.y(2:end); Manif.points.z(2:end)]);
rem_deltamin.idx=find(rem_deltamin.delta<manif.growinf.deltamin); %candidates to be removed (idx)         
rem_deltamin.rem=[]; %points we are going to remove

%we want to see if the idx are consecutive or not
idx1=rem_deltamin.idx(1:end-1);
idx2=rem_deltamin.idx(2:end);
diff=idx2-idx1;

[~,breaks]=find(diff~=1); %find consecutive indexes. If all are consecutive breaks=0

if numel(breaks)==0 %if all indexes are consecutive, start from the first index beggining
    breaks=rem_deltamin.idx(1);
end

for i=1:numel(breaks)
    if breaks(1:end)==rem_deltamin.idx(1) %if there are no breaks, take all indexes
        part=rem_deltamin.idx;
    elseif i==1
        part=[rem_deltamin.idx(1:breaks(i)) rem_deltamin.idx(breaks(i))+1];
    else
        part=[rem_deltamin.idx(breaks(i-1)+1:breaks(i)) rem_deltamin.idx(breaks(i))+1];
    end

    t=cumsum(rem_deltamin.delta(part)); %To know the distances of each point from the first one
    tt=0:manif.growinf.deltamin:t(end); % this is the (smallest)ideal distance of points

    % I keep the nearest points to the ideal from the mesh
    % Here 't' is the distance, 'part' is the corresponding indices of those distances and 'tt' is the ideal.
    keep=interp1(t,part,tt(2:end),'nearest'); 
    %----------------

    % we always keep points at infinity in z coordinate, for avoiding trimming the corners
    inf=find(abs(Manif.points.z(part))>0.999);
    [~,maxx]=max(abs(Manif.points.z(part)));
    keep=unique([keep part(intersect(inf,maxx))]);

    rem_deltamin.rem=[ rem_deltamin.rem part(~ismember(part,keep))]; %points we are going to remove, exept super near infinity


end

Manif.growinf.runinf.rem_deltamin=Manif.growinf.runinf.rem_deltamin+numel(rem_deltamin.rem);
Manif.points.x(rem_deltamin.rem)=[]; 
Manif.points.y(rem_deltamin.rem)=[]; 
Manif.points.z(rem_deltamin.rem)=[]; 
    
%-------------- END Removing points depending on Deltamin

Manif.points.arc = arclength(Manif.points);

fprintf('\n arclength %.0f\n',Manif.points.arc(end))
fprintf('\n');   
fprintf('\n');

Manif.growinf.runinf.npoints_initial_final=[Manif.growinf.runinf.npoints_initial_final numel(Manif.points.x)];
Manif.growinf.runinf.arc_initial_final=[Manif.growinf.runinf.arc_initial_final Manif.points.arc(end)];
Manif.growinf.runinf.time=toc;




fprintf('\n Elapsed time is %.3f seconds\n\n',Manif.growinf.runinf.time)
fprintf('\n %i initial points, arclength %.0f\n',Manif.growinf.runinf.npoints_initial_final(1),Manif.growinf.runinf.arc_initial_final(1))
fprintf(' %i final points, arclength %.0f \n',numel(Manif.points.x),Manif.points.arc(end)) %76800 longer
fprintf('   * %i points removed \n',Manif.growinf.runinf.rem_deltamin+Manif.growinf.runinf.rem_nan+Manif.growinf.runinf.rem_inf) 
fprintf('   * %i points added from deltamax \n',Manif.growinf.runinf.add_deltamax) 
fprintf('   * %i points added from alpha \n',Manif.growinf.runinf.add_alphamax) 
fprintf('   * %i points added from delta*alpha \n',Manif.growinf.runinf.add_deltalphamax) 
%% FUNCTIONS

function interp = makima3D(points,t,tt)
% get interpolation points
% t: parametrization of the points
% tt: parametrization of the interpolated points
interp=struct();

interp.x = interp1(t,points.x,tt,'makima','extrap');
interp.y = interp1(t,points.y,tt,'makima','extrap');
interp.z = interp1(t,points.z,tt,'makima','extrap');
end 

%----------------

function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)
arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength

%----------------

function alpha=angles(p0,p1,p2)
%angle between p0p1 and p1p2
n1 = (p1 - p2) / norm(p1 - p2);  % Normalized vectors
n2 = (p0 - p1) / norm(p0 - p1);
alpha = atan2(norm(cross(n1, n2)), dot(n1, n2)); %gives value from 0 to pi
end

%----------------

function Anew=insert(A,B,ind)
% Anew: new vector with the new values
% A: Old vector
% B: vector with new values
%ind: index to insert values after this row


    % Preallocate output
    Anew = zeros(1,numel(A)+numel(B));

    % Find indices for old data
    addRows = ismember(1:numel(A), ind);
    oldDataInd = (1:numel(A)) + cumsum([0, addRows(1:end-1)]);

    % Add in old data
    Anew(oldDataInd) = A;

    % Find indices for new data
    newDataInd = (1:length(ind)) + ind;

    % Add in new data
    Anew(newDataInd) = B;
end

end





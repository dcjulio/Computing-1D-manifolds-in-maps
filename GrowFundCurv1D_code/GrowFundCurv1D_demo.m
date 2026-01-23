%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions_v3.0.0');

%% Options

    %--- Information of the system
    opts.thesystem=StdHenon3D; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', -0.3, 'xi', 0.8); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 2000; % What is the approximate arclength of the manifold
    
    %--- Number of iterations used to compute the manifold
    opts.max_funditer = 100; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7; % pos value is positive branch and viceversa


%% Computing the manifold: Ws(pmin) orientation-preserving
    opts.name_fixpoint='pplu'; % name of the fixed point associated to the manifold (has to match with the names defined in StdHenon3D)
    opts.stability='Umanifold';

    opts.branch='pos'; %which branch (Optional. Could be done with the sign of the init_step instead)
    manif=GrowFundCurv1D(opts);

    %% Computing the other branch (if orientation preserving)
    % manif=add_branch(manif, opts, 'neg');

%% Computing intersection points
    angle=pi/4; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif=inter_plane(manif,angle);

%% Plot
    manifplot(manif);

%% Epsilon pseudo orbit (orientation-preserving)

branch    = opts.branch; %what branch is intersecting the plane

if isfield(manif,'inter') && numel(manif.inter.points.(branch).idx)>0
    
    idxpoint  = manif.inter.points.(branch).idx(end); % for example, the last intersection point
    orbit     = eps_pseudo_orbit(manif, idxpoint, branch); % computing the pseudo orbit of the intersection point

% Plotting the eps-pseudo-orbit
     hold on
     plot3(orbit.x([1,end]),orbit.y([1,end]),orbit.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
     plot3(orbit.x,orbit.y,orbit.z,'ko--','LineWidth',1.2) %epsilon orbit
end



%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    %--- Information of the system
    opts.thesystem=StdHenon3D; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength=2000; % What is the approximate arclength of the manifold
    
    %--- Number of iterations used to compute the manifold
    opts.max_funditer=60; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7; % pos value is positive branch and viceversa


%% Computing the manifold: Ws(pmin) orientation-preserving
    opts.name_fixpoint='pmin'; % name of the fixed point associated to the manifold (has to match with the names defined in StdHenon3D)
    opts.stability='Smanifold';

    opts.branch='pos'; %which branch (Optional. Could be done with the sign of the init_step instead)
    manif=GrowFundCurv1D(opts);

    %% Computing the other branch
    manif=add_branch(manif, opts, 'neg');

%% Computing intersection points
    angle=-3*pi/4; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif=inter_plane(manif,angle);

%% Plot
    manifplot(manif);

%% Adding the path
clear all
addpath('./GrowFundCurv1D_functions');

%% Options
%--- Loading and saving the initial segment of manifold
load('pmin_a4.2_b0.3_xi1.20.mat')
opts.init_segment=struct('x',data_x,'y',data_y,'z',data_z);

%--- Information of the system
opts.thesystem=StdHenon3D; %what is the name of the class with the info of the diffeo
opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); %the parameter values
opts.name_fixpoint='pmin'; % fixed point associated to the manifold
opts.branch='pos'; %what side of the manifold I am computing. It is just a name. Examples: 'neg', 'pos', ''. 

%--- Number of iterations used to compute the manifold
opts.funditer=6; %whow many times the algorithm iterates the fundamental domain (multiple of mapiter).
opts.mapiter=2; % how many iterations of the map are used (if non-orientable has to be even)

%--- Accuracy parameters (default)
%opts.accpar.alphamax=0.3;
%opts.accpar.deltalphamax=0.001; 
%opts.accpar.deltamin=0.000001;
%opts.accpar.deltamax=0.01;    

%% Initializing the manifold structure
manif = init_manif(opts);

%% Computing the manifold
manif=GrowFundCurv1D(manif,opts);

%% Computing intersection points
opts.angle=-3*pi/4; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
manif=inter_plane(manif,opts);

%% Plotting the data
f1=manifplot(manif);

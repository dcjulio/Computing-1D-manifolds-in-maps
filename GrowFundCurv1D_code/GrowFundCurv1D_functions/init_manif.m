function manif = init_manif(opts)

%---- manif.name: Name of the manifold---% 
% for example Ws_pmin_pos is the stable manifolds of the fixed point pmin, the branch to the positive values
%
%---- manif.orientability: Orientability of the manifold---% 
%
%---- manif.fixp: Information of the fixed points associated to the manifold---% 
%
%---- manif.stab: Stability of the manifold---% 
%
%---- manif.points: Coordinates of the manifold ---% 
%
%---- manif.inf_sys: contains general info about the map ---% 
% inf_sys.par: The parameter values
% inf_sys.fixp: the fixed points with their  eigensystem, orientability, stability, etc
%
%---- manif.growinf: information for the algorithm---% 
%
%%
%the map function where the system is defined (example StdHenon3D)
thesystem=opts.thesystem;   

%% Initializate field names and the structure 'manif'

%names of the fields
names = {
    'name'
    'orientability'
    'fixp'   % fixpoint
    'stab'   % stability
    'points'
    'inf_sys'
    'growinf'    % options for computation
    }; 

manif=struct();
n=numel(names);

for k=1:n
    manif.(names{k})=[];
end

%% General info of the map
manif.inf_sys.par=opts.par;  % parameters
manif.inf_sys.fixp=thesystem.fixpoints(opts); % fixed points


% eigensystem, stability and orientability of all the fixed points
fixp_names=fieldnames(manif.inf_sys.fixp);
fixpinfo=manif.inf_sys.fixp;
for k=1:numel(fixp_names)

    % eigensystem
    [fixpinfo.(fixp_names{k}).eigval,fixpinfo.(fixp_names{k}).eigvec]=eigensystem(fixpinfo.(fixp_names{k}),opts);

    % recognize the stability of the one dimensional manifold
    fixpinfo.(fixp_names{k}).stab=stability(fixpinfo.(fixp_names{k}).eigval);
    
    % orientability of the one dimensional manifold
    [fixpinfo.(fixp_names{k}).orientability]=orientability(fixpinfo.(fixp_names{k}).stab,fixpinfo.(fixp_names{k}).eigval);
end
%saving the info
manif.inf_sys.fixp=fixpinfo;
%% Name of the manifold. Example: Ws_pmin_pos

name_fixpoint=opts.name_fixpoint; 
branch=opts.branch;

% defining the name of the manifold
if strcmp(branch,'')
    manif.name = sprintf('W%s_%s', lower(manif.inf_sys.fixp.(name_fixpoint).stab(1)),name_fixpoint);
else
    manif.name = sprintf('W%s_%s_%s', lower(manif.inf_sys.fixp.(name_fixpoint).stab(1)),name_fixpoint,branch);
end
  
%% Orientability of the manif
manif.orientability=manif.inf_sys.fixp.(name_fixpoint).orientability; 	% orientability of the manifold

%% Coordinate and eigensystem of the fixed point associated to the manifold
manif.fixp.(name_fixpoint).x=manif.inf_sys.fixp.(name_fixpoint).x; 	% values of fixpoint coord x
manif.fixp.(name_fixpoint).y=manif.inf_sys.fixp.(name_fixpoint).y; 	% values of fixpoint coord y
manif.fixp.(name_fixpoint).z=manif.inf_sys.fixp.(name_fixpoint).z; 	% values of fixpoint coord z
    
manif.fixp.eigsys.eigval=manif.inf_sys.fixp.(name_fixpoint).eigval;    % eigenvalue of fixpoint
manif.fixp.eigsys.eigvec=manif.inf_sys.fixp.(name_fixpoint).eigvec;    % eigenvector of fixpoint

%% Stability of the manifold
manif.stab=manif.inf_sys.fixp.(name_fixpoint).stab;               % stability

%% Initial segment of the manifold
manif.points=struct('x',opts.init_segment.x,'y',opts.init_segment.y,'z',opts.init_segment.z);
manif.points.arc=arclength(manif.points);

%% algorithm information
manif.growinf.mapiter=opts.mapiter; % number of iteration of the map
manif.growinf.funditer=opts.funditer; % number of iteration of the algorithm
%% Accuracty conditions

% default acc conditions
manif.growinf.alphamax=0.3;
manif.growinf.deltalphamax=0.001; 
manif.growinf.deltamin=0.000001;
manif.growinf.deltamax=0.01;    


%rewrite them if we other values are defined
if isfield(opts,'accpar')
    if isfield(opts.accpar,'alphamax')
        manif.growinf.alphamax=opts.accpar.alphamax;
    end
    if isfield(opts.accpar,'deltalphamax')
        manif.growinf.deltalphamax=opts.accpar.deltalphamax;
    end
    if isfield(opts.accpar,'deltamin')
        manif.growinf.deltamin=opts.accpar.deltamin;
    end
    if isfield(opts.accpar,'deltamax')
        manif.growinf.deltamax=opts.accpar.deltamax;
    end
end

    

%----------------------------------------------
%-------------- FUNCTIONS ---------------------
%----------------------------------------------

% > -------- eigenvalues and eigenvectors
function [eigval, eigvec]=eigensystem(fixpoint,opts)

    
    syms x y z
    points.x=x;
    points.y=y;
    points.z=z;
    
    map=StdHenon3D.ff(points,opts);
    
    J=jacobian([map.x, map.y, map.z], [x, y, z]);
    
    %decompactify fixed points
    fixpoint=StdHenon3D.decompactify(fixpoint);
    
    x=fixpoint.x;
    y=fixpoint.y;
    z=fixpoint.z;

    [eigvec,D]=eig(J);
    eigval=diag(D);
    
    eigvec=double(subs(eigvec));
    eigval=double(subs(eigval));
    
end

%----------------------------------------------
%----------------------------------------------
%----------------------------------------------

% > -------- stability
function Stab=stability(eigval)
    
    %if there is only one eigenvalue bigger than 1, then the one
    %dimensional manifold is Unstable
    if sum(abs(eigval)>1)==1
        Stab='Umanifold';  % 'Umanifold', 'Smanifold'
    else
        Stab='Smanifold';
    end
    
end
%----------------------------------------------
%----------------------------------------------
%----------------------------------------------

% > -------- orientability
function [orientability]=orientability(Stab,eigval)

    % check the sign of the eigenvalue
    if strcmp(Stab,'Umanifold')
        if eigval(abs(eigval)>1)>0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end
    
    if strcmp(Stab,'Smanifold')
        if eigval(abs(eigval)<1)>0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end
    
end  

% > -------- oarclength
function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)

arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength
    
end
function manif = init_manif(opts)

%---- manif.name: Name of the manifold---% 
% for example Ws_pmin_pos is the stable manifolds of the fixed point pmin, the branch to the positive values
%
%---- manif.orientability: Orientability of the manifold---% 
%
%---- manif.fixp: Information of the fixed points associated to the manifold---% 
%
%---- manif.stability: Stability of the manifold---% 
%
%---- manif.points: Coordinates of the manifold ---% 
%
%---- manif.system_info: contains general info about the map ---% 
% system_info.par: The parameter values
% system_info.fixp: the fixed points with their  eigensystem, orientability, stability, etc
%
%---- manif.grow_info: information for the algorithm---% 
%
%%
%the map function where the system is defined (example StdHenon3D)
thesystem=opts.thesystem;   
%% Initializate field names and the structure 'manif'

%names of the fields
names = {
    'name'
    'orientability'
    'stability'   % stability
    'dimension'
    'fixp'   % fixpoint
    'system_info'
    'grow_info'    % options for computation
    }; 

manif=struct();
n=numel(names);

for k=1:n
    manif.(names{k})=[];
end

%% Default accuracty conditions

% default acc conditions
manif.grow_info.alphamax=0.3;
manif.grow_info.deltalphamax=10^(-3); 
manif.grow_info.deltamin=10^(-6);
manif.grow_info.deltamax=10^(-2);    
manif.grow_info.init_step=10^(-7);  

%% Update accuracy conditions if needed
manif.grow_info.thesystem = thesystem;

%rewrite them if we other values are defined
if isfield(opts,'accpar')
    if isfield(opts.accpar,'alphamax')
        manif.grow_info.alphamax=opts.accpar.alphamax;
    end
    if isfield(opts.accpar,'deltalphamax')
        manif.grow_info.deltalphamax=opts.accpar.deltalphamax;
    end
    if isfield(opts.accpar,'deltamin')
        manif.grow_info.deltamin=opts.accpar.deltamin;
    end
    if isfield(opts.accpar,'deltamax')
        manif.grow_info.deltamax=opts.accpar.deltamax;
    end
    if isfield(opts.accpar,'init_step')
        manif.grow_info.init_step=opts.accpar.init_step;
    end
end


%% General info of the map
manif.system_info.par=opts.par;  % parameters
manif.system_info.fixp=thesystem.fixpoints(opts); % fixed points. 'thesystem' store info about the map

% eigensystem, stability and orientability of all the fixed points
fixp_names=fieldnames(manif.system_info.fixp);
fixpinfo=manif.system_info.fixp;

for k=1:numel(fixp_names)

    % eigensystem
    [fixpinfo.(fixp_names{k}).eigval,fixpinfo.(fixp_names{k}).eigvec]=eigensystem(fixpinfo.(fixp_names{k}),opts);
    eigval=fixpinfo.(fixp_names{k}).eigval;
    eigvec=fixpinfo.(fixp_names{k}).eigvec;
    fixpinfo.(fixp_names{k}).eigvec=fixpinfo.(fixp_names{k}).eigvec;


    % computes the dimension and orientation properties of the stable manifold
    fixpinfo.(fixp_names{k}).Smanifold.dimension=sum(abs(eigval)<1);
    fixpinfo.(fixp_names{k}).Smanifold.orientability=orientability(eigval,'Smanifold');
    fixpinfo.(fixp_names{k}).Smanifold.eigval=eigval(abs(eigval)<1);
    fixpinfo.(fixp_names{k}).Smanifold.eigvec=eigvec(abs(eigval)<1,:);


    % computes the dimension and orientation properties of the unstable manifold
    fixpinfo.(fixp_names{k}).Umanifold.dimension=sum(abs(eigval)>1);
    fixpinfo.(fixp_names{k}).Umanifold.orientability=orientability(eigval,'Umanifold');
    fixpinfo.(fixp_names{k}).Umanifold.eigval=eigval(abs(eigval)>1);
    fixpinfo.(fixp_names{k}).Umanifold.eigvec=eigvec(abs(eigval)>1,:);
    
end
%saving the info
manif.system_info.fixp=fixpinfo;

 
%% Stability and orientability and dimension of the manif to compute

name_fixpoint=opts.name_fixpoint; 

manif.stability=opts.stability;
manif.orientability=manif.system_info.fixp.(name_fixpoint).(opts.stability).orientability; 	% orientability of the manifold
manif.dimension=manif.system_info.fixp.(name_fixpoint).(opts.stability).dimension;

%if the field branch is defined, then follow that definition to know which branch to compute
if isfield(opts,'branch')
   if strcmp(opts.branch,'pos')
       manif.grow_info.init_step=abs(manif.grow_info.init_step);
   elseif strcmp(opts.branch,'neg')
       manif.grow_info.init_step=-abs(manif.grow_info.init_step);
   end
end

%% Name of the manifold. Example: Ws_pmin_pos

% defining the name of the manifold
manif.name = sprintf('W%s_%s', lower(manif.stability(1)),name_fixpoint);  

%% Coordinate and eigensystem of the fixed point associated to the manifold
manif.fixp=manif.system_info.fixp.(name_fixpoint); 
manif.fixp.name=name_fixpoint;

%% algorithm information
manif.grow_info.stability=manif.stability; % stability of the manifold
manif.grow_info.orientability=manif.orientability; % orientability of the manifold
manif.grow_info.dimension=manif.dimension; % dimension of the manifold


manif.grow_info.eigval=manif.fixp.(manif.stability).eigval; 
manif.grow_info.eigvec=manif.fixp.(manif.stability).eigvec; 

manif.grow_info.max_funditer=opts.max_funditer; % number of iteration of the algorithm
manif.grow_info.max_funditer=opts.max_funditer; % number of iteration of the algorithm
manif.grow_info.user_arclength=opts.user_arclength;

    

%----------------------------------------------
%-------------- FUNCTIONS ---------------------
%----------------------------------------------

% > --------  eigenvalue and eigenvector in compactified coordinates
function [eigval, eigvec]=eigensystem(fixpoint,opts)

    syms x y z
    points=struct('x',x,'y',y,'z',z);
   
    system=opts.thesystem;

    %Jacobian J_f of the original system (uncompactified)
    F=system.ff(points,opts);
    JF=jacobian([F.x, F.y, F.z], [x, y, z]);
    fixp=system.decompactify(fixpoint); %original fixed point

    %Jacobian J_t of the compactification
    T=system.compactify(points);
    JT=jacobian([T.x, T.y, T.z], [x, y, z]);
    

    %evaluate the jacobians in the fixed points
    x=fixp.x;
    y=fixp.y;
    z=fixp.z;
    JFp=double(subs(JF));

    JTp=double(subs(JT));

    %Compute the eigenval and eigenvec of the original system JF(p)
    [eigvecF,D]=eig(JFp);
    eigval=diag(D)';
    
    %Transform the eigenvec to the xompaxtified coordinates
    eigvec=JTp*eigvecF;
    eigvec=normc(eigvec); %normalize each column

    eigvec=eigvec'; %each row is a vector
    eigval=eigval';
end

%----------------------------------------------
%----------------------------------------------
%----------------------------------------------


% > -------- orientability
function [orientability]=orientability(eigval,Stab)

    if strcmp(Stab,'Smanifold')
        if prod(eigval(abs(eigval) < 1))>0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end


    if strcmp(Stab,'Umanifold')
        if prod(eigval(abs(eigval) > 1))>0
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
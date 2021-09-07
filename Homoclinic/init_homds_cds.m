function init_homds_cds(odefile,bt, ap, options)
%% Auxiliary function to initialize global variables
%  homds and cds for the continuation of the homclinic solution.
%  These global variables are used throughout MatCont

ntst = options.ntst;
ncol = options.ncol;
extravec = options.extravec;

global homds
homds = [];

global cds

if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end

cds.curve = @homoclinic;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};

% set parameters for Newton corrections in newtcorr.m if not set yet
if isempty(cds.options.MaxCorrIters)
    cds.options.MaxCorrIters = 10;
end
if isempty(cds.options.MaxNewtonIters)
    cds.options.MaxNewtonIters = 3;
end
if isempty(cds.options.MoorePenrose)
    cds.options.MoorePenrose = 1;
end
if isempty(cds.options.VarTolerance)
    cds.options.VarTolerance = 1e-08;
end
if isempty(cds.options.FunTolerance)
    cds.options.FunTolerance = 1e-08;
end

if isempty(cds.options.Increment)
    cds.options.Increment=1.0e-05;
end

homds.odefile = odefile;
func_handles = feval(homds.odefile);
if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
else symord = [];
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
else symordp = [];
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
%cds.symjac = 1;
%cds.symhess = 0;

homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
%
siz = size(func_handles,2);

if siz > 9
    j=1;
    for k=10:siz
        homds.user{j}= func_handles{k};
        j=j+1;
    end
else homds.user=[];
end

homds.nphase = length(bt.x);
homds.ActiveParams = ap;
homds.P0 = bt.par;
homds.extravec = extravec;

Hom_set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);

homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.pars = homds.ncoords+(1:3);
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;
homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.bialt_M1 = [];
homds.bialt_M2 = [];
homds.bialt_M3 = [];
homds.bialt_M4 = [];
homds.multipliers = nan;
homds.monodromy = [];
homds.multi_r1 = [];
homds.multi_r2 = [];
homds.ups = [];
homds.vps = [];
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;

homds.HTPstep = 0;

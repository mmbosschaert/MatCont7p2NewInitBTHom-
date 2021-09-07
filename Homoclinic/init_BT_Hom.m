function [x1, v1, x1_pred, v1_pred] = init_BT_Hom(odefile, bt, ap, options)

%% Check input
if length(ap) ~= 2
    error('2 free system parameters are needed');
elseif options.extravec==zeros(1,3)
    error('at least 1 free homoclinic parameter is needed');
elseif options.extravec==ones(1,3)
    error('at most 2 free homoclinic parameters are needed');
end

%% Initialize cds and homds
init_homds_cds(odefile, bt, ap, options);
global homds cds

for i=1:options.MaxNumCorrections
    %% Select method
    switch options.method
        case 'orbital'
            [ups, p, dalphadepsilon] = init_BT_Hom_orbital(odefile, bt, ap, options);
        case 'orbitalv2'
            [ups, p, dalphadepsilon] = init_BT_Hom_orbital_v2(odefile, bt, ap, options);
        case 'LP'
            [ups, p, dalphadepsilon] = init_BT_Hom_LP(odefile, bt, ap, options);
        case 'LPHypernormalForm'
            [ups, p, dalphadepsilon] = init_BT_Hom_without_e_b1(odefile, bt, ap, options);
        case 'RegularPerturbation'
            [ups, p, dalphadepsilon] = init_BT_Hom_RP(odefile, bt, ap, options);
        case 'RegularPerturbationL2'
            [ups, p, dalphadepsilon] = init_BT_Hom_RP_L2(odefile, bt, ap, options);
        otherwise
            error('init_BT_Hom: No method selected');
    end

    %% The initial eps0 & eps1 :
    homds.eps0 = norm(ups(:,  1) - homds.x0);
    homds.eps1 = norm(ups(:,end) - homds.x0);
    if options.messages
        fprintf('The initial distance eps0: %g\n', homds.eps0);
        fprintf('The initial distance eps1: %g\n', homds.eps1);
    end

    %% Dimension of stable & unstable subspace
    Hom_calc_weights;
    Jac = cjac(homds.func,homds.Jacobian,homds.x0,num2cell(p),homds.ActiveParams);
    D = eig(Jac);
    homds.nneg  = sum(real(D) < 0);
    homds.npos  = homds.nphase-homds.nneg;
    homds.Ysize = homds.nneg*homds.npos;

    %% Compose x1
    % i. Orbit
    x1 = reshape(ups,[],1);
    v = [];
    [x1,~]=Hom_new_mesh(x1,v,options.ntst,options.ncol);
    % ii. EQUILIBRIUM COORDINATES
    x1 = [x1; homds.x0];
    % iii. TWO FREE PARAMETERS
    x1 = [x1; homds.P0(homds.ActiveParams)];
    % iv. EXTRA FREE PARAMETER (FREE HOMOCLINIC PARAMETER)
    extravec = [homds.T; homds.eps0; homds.eps1];
    x1 = [x1; extravec(find(homds.extravec))];
    for i=1:homds.nneg
        x1 = [x1; zeros(homds.npos,1)];
    end
    for i=1:homds.npos
        x1 = [x1; zeros(homds.nneg,1)];
    end

    %% Assign values to homoclinic fields
    % a) YS AND YU, INITIALIZED TO 0
    %
    homds.YS = zeros(homds.npos,homds.nneg);
    homds.YU = zeros(homds.nneg,homds.npos);
    % b) THIRD PARAMETER = UNSTABLE_FLAG:-
    % i) 1 IF WE WANT THE UNSTABLE SPACE.
    % ii) 0 IF WE WANT THE STABLE ONE.
    [QU, ~] = computeBase(Jac,0,homds.npos);
    [QS, ~] = computeBase(Jac,1,homds.nneg);
    homds.oldStableQ = QS;
    homds.oldUnstableQ = QU;
    homds.ups = [];
    homds.ndim = length(x1);
    cds.ndim = homds.ndim; % here was a major typo cd.ndim=homds.ndim

    %% Initial tangent vector
    ups0 = reshape(x1(1:size(ups,1)*size(ups,2),:),homds.nphase,homds.tps);
    pp1 = num2cell(p);
    for i=1:homds.tps
        homds.upoldp(:,i) = 2*homds.T*feval(homds.func, 0, ups0(:,i), pp1{:});
    end
    BTHomjac = BVP_Hom_jac(homds.func,x1(1:size(ups,1)*size(ups,2),1),...
        homds.x0,p,homds.T,homds.eps0,homds.eps1,homds.YS,homds.YU);
    [Q,~] = qr(full(BTHomjac)');
    v1 = full(Q(:,end));
    % v2 = null(full(BTHomjac));

    %% Change, if necessary, the direction to the tangent vector
    % we should start in the direction away from the Bogdanov-Takens point
    if dalphadepsilon(1)*v1(homds.ncoords+ homds.nphase + 1) < 0
       v1 = -v1;
    end

    x1_pred = x1;
    v1_pred = v1;
    %% Correct with Newton
    if options.correct
        [x1, v1] = newtcorr(x1, v1);
        if isempty(x1)
            % try with smaller perturbationparameter if newton didn't converge
           options.perturbationparameter = options.perturbationparameter/2;
        else
            break;
        end
    else
        break;
    end
end

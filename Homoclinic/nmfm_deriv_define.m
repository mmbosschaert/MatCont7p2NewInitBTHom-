function F = nmfm_deriv_define(odefile, bt, ap)

global homds;
global cds;

x = bt.x;
p = bt.par;

func_handles = feval(odefile);
symord = 0;
symordp = 0;
%
if     ~isempty(func_handles{9}),   symord = 5;
elseif ~isempty(func_handles{8}),   symord = 4;
elseif ~isempty(func_handles{7}),   symord = 3;
elseif ~isempty(func_handles{5}),   symord = 2;
elseif ~isempty(func_handles{3}),   symord = 1;
end
if     ~isempty(func_handles{6}),   symordp = 2;
elseif ~isempty(func_handles{4}),   symordp = 1;
end


if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end

cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;
%
homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};
%
cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];
cds.ndim = length(xp);
%
homds.x0 = bt.x;
%
pcell = num2cell(p);

hessIncrement=(cds.options.Increment)^(3.0/4.0);
ten3Increment=(cds.options.Increment)^(3.0/5.0);

if (cds.options.SymDerivative>=3)
  hess= chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,ap);
  tens=ctens3(homds.func,homds.Jacobian,homds.Hessians,homds.Der3,homds.x0,pcell,ap);
else
  hess = [];
  tens = [];
end
F.A  =  cjac(homds.func,homds.Jacobian, homds.x0,pcell, ap);
F.J1 = cjacp(homds.func,homds.JacobianP,homds.x0,pcell,ap);
F.B  = @(v1,v2) multilinear2(homds.func,hess,v1,v2,homds.x0,pcell,hessIncrement);
F.C  = @(v1,v2,v3) multilinear3(homds.func,tens,v1,v2,v3,homds.x0,pcell,ten3Increment);

%% define multilinarforms file is exists
multilinarforms_file= [char(odefile), '_multilinearforms'];
if isOctave
    multilinarforms_file = regexprep(multilinarforms_file, '@', '');
end
multilinarforms = [];
if exist(multilinarforms_file, 'file')
    multilinarforms = eval(multilinarforms_file);
end

par0 = p*0; 
dp = @(vals) matcont_array_insert(par0, ap, vals);
if isfield(multilinarforms,'J2')
    F.J2 = @(p1, p2) multilinarforms.J2(x', p', dp(p1), dp(p2));
else
    % J2FD
    for i=ap
        pa1 = pcell; pa1{i} = pa1{i}-cds.options.Increment;
        pa2 = pcell; pa2{i} = pa2{i}+cds.options.Increment;
        Hjp2 = cjacp(homds.func,homds.JacobianP,homds.x0,pa2,ap);
        Hjp1 = cjacp(homds.func,homds.JacobianP,homds.x0,pa1,ap);
        tmpJ2FD(:,:,i) = Hjp2 - Hjp1;
    end
    tmpJ2FD = tmpJ2FD/(2*cds.options.Increment);
    J2FD = tmpJ2FD;
    if size(J2FD,3) > length(ap)
        J2FD = J2FD(:,:,ap);
    end
    F.J2 = @(p1,p2) tensor2op(J2FD,p1,p2,2);
end
if isfield(multilinarforms,'B1')
    F.B1 = @(v1, v2, p1) multilinarforms.B1(x', p', v1, v2, dp(p1));
else
    % B1FD
    for i=ap
        pa1= pcell; pa1{i} = pa1{i}-cds.options.Increment;
        pa2= pcell; pa2{i} = pa2{i}+cds.options.Increment;
        Hp2=chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pa2,ap);
        Hp1=chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pa1,ap);
        %
        B1FD(:,:,:,i)=Hp2-Hp1;
    end
    B1FD=B1FD/(2*cds.options.Increment);
    if size(B1FD,4) > length(ap)
        B1FD=B1FD(:,:,:,ap);
    end
    F.B1 = @(v1,v2,p1) [tensor2op(B1FD(:,:,:,1),v1,v2,homds.nphase),...
        tensor2op(B1FD(:,:,:,2),v1,v2,homds.nphase)]*p1;
end
if isfield(multilinarforms,'A1')
    F.A1 = @(v1, p1) multilinarforms.A1(x', p', v1, dp(p1));
else
    A1FD   = chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pcell,ap);
    if size(A1FD,3) > length(ap)
        A1FD = A1FD(:,:,ap);
    end
    F.A1 = @(v1,p1) tensor2op(A1FD,v1,p1,2);
end
if isfield(multilinarforms,'A2')
    F.A2 = @(v1, p1, p2) multilinarforms.A2(x', p', v1, dp(p1), dp(p2));
else
    % A2FD
    for i=ap
        pa1= pcell; pa1{i} = pa1{i}-cds.options.Increment;
        pa2= pcell; pa2{i} = pa2{i}+cds.options.Increment;
        Hp2=chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pa2,ap);
        Hp1=chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pa1,ap);
        A2FD(:,:,:,i)=Hp2-Hp1;
    end
    A2FD=A2FD/(2*cds.options.Increment);
    if size(A2FD,4) > length(ap)
        A2FD=A2FD(:,:,:,ap);
    end
    F.A2 = @(v1,p1,p2) [tensor2op(A2FD(:,:,:,1),v1,p1,2), tensor2op(A2FD(:,:,:,2),v1,p1,2)]*p2;
end
if isfield(multilinarforms,'J3')
    F.J3 = @(p1, p2, p3) multilinarforms.J3(x', p', dp(p1), dp(p2), dp(p3));
else
    F.J3 = @(p1, p2, p3) D0D3F(homds.func, p1, p2, p3, homds.x0, p, ten3Increment,ap);
end

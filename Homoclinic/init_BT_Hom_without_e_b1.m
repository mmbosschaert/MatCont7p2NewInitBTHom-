function [ups, p, dalphadepsilon] = init_BT_Hom_without_e_b1(odefile, bt, ap, options)

global homds cds

%% 3. Compute normal form coefficients
bt = BT_nmfm_without_e_b1(odefile, bt, ap, options);
a  = bt.nmfm.a;
b  = bt.nmfm.b;
d  = bt.nmfm.d;
e  = bt.nmfm.e;
a1 = bt.nmfm.a1;
b1 = bt.nmfm.b1;
K10   = bt.nmfm.K10;
K01   = bt.nmfm.K01;
K02   = bt.nmfm.K02;
K11   = bt.nmfm.K11;
K03   = bt.nmfm.K03;
q0    = bt.nmfm.q0;
q1    = bt.nmfm.q1;
H0010 = bt.nmfm.H0010;
H0001 = bt.nmfm.H0001;
H2000 = bt.nmfm.H2000;
H1100 = bt.nmfm.H1100;
H0200 = bt.nmfm.H0200;
H1010 = bt.nmfm.H1010;
H1001 = bt.nmfm.H1001;
H0110 = bt.nmfm.H0110;
H0101 = bt.nmfm.H0101;
H0002 = bt.nmfm.H0002;
H0011 = bt.nmfm.H0011;
H3000 = bt.nmfm.H3000;
H2100 = bt.nmfm.H2100;
H1101 = bt.nmfm.H1101;
H2001 = bt.nmfm.H2001;
H0003 = bt.nmfm.H0003;
H1002 = bt.nmfm.H1002;
H0102 = bt.nmfm.H0102;
if options.messages
    disp('BT normal form coefficients:')
    fprintf('a=%d,\t b=%d\n', a, b);
end

%% Set amplitude and TTolerance and HigherOrderTimeReparametrization
if ~isempty(options.amplitude)
    amplitude = options.amplitude;
    eps  = sqrt(amplitude*abs(a)/6);
else
    eps = a/b*options.perturbationparameter;
    amplitude = eps^2/abs(a)*6;
end
if ~isempty(options.TTolerance)
    TTolerance = options.TTolerance;
else
    TTolerance = amplitude*options.amplitudeTToleranceRatio;
end
if ~isempty(options.HigherOrderTimeReparametrization)
    HigherOrderTimeReparametrization = options.HigherOrderTimeReparametrization;
else
    HigherOrderTimeReparametrization = 1;
end

%% 8. Initial cycle
% a) The initical approximation of the parameters
beta1 = -4/a*eps^4;
tau0 = 10/7;
tau2 = 1/a*(100/49*b1 - 4*e/b) + 1/a^2*(288/2401*b^2-50*b*a1/49 + 146/49*d);
beta2 = b/a*(tau0 + tau2*eps^2)*eps^2;
alpha= K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2 + K11*beta1.*beta2 + ...
    1/6*K03*beta2^3;
p = bt.par;
p(ap) = p(ap) + alpha;  % tmpfreep=tmpfreep+alpha;
homds.P0 = p;

%% b) THE INITIAL HALF-RETURN TIME VALUE T
% |w0(-inf,eps)-w0(T,eps)|=TTolerance=>(6*norm(q0)*eps^2/|a|)*tanh(eps*T)^2=k
% where norm(q0)=1.
% notice that TTolerance is Tolerance (and k) in Al-Hdabat et.al. 2016.
% notice furthermore that T was calculated incorrect in the code
homds.T = 1/eps*asech(sqrt(TTolerance/amplitude));
if options.messages
    fprintf('The initial perturbation parameter:  %d\n', eps);
    fprintf('The initial amplitude: %g\n', amplitude);
    fprintf('T: %g\n', homds.T);
end

%% c) THE INITIAL APPROXIMATION OF CYCLE
t = (2*homds.finemsh - 1) * (homds.T);
% transformation of time \xi to s
xi1 = @(s) -6 .*b .*log(cosh(s))/7/a;
xi2 = @(s) (70 .*a1 .*b .*s-72 .*b.^2 .*s+196 .*d .*s+90 .*b.^2 .*tanh(s)-441 .*d .*tanh(s)+144 .*b.^2 .*log(cosh(s)) .*tanh(s))/(196 .*a.^2);
xi3 = @(s) (1/(4802 .*a.^3)) .*(4410 .*a1 .*b.^2 .*log(cosh(s))-522 .*b.^3 .*log(cosh(s))-5880 .*a .*b .*b1 .*log(cosh(s))+12201 .*b .*d .*log(cosh(s))-1512 .*b.^3 .*log(cosh(s)).^2-252 .*b.^3 .*sech(s).^2+10290 .*b .*d .*sech(s).^2-9604 .*a .*e .*sech(s).^2-1470 .*a1 .*b.^2 .*s .*tanh(s)+1512 .*b.^3 .*s .*tanh(s)-4116 .*b .*d .*s .*tanh(s)-1890 .*b.^3 .*tanh(s).^2+9261 .*b .*d .*tanh(s).^2-1134 .*b.^3 .*log(cosh(s)) .*tanh(s).^2-9261 .*b .*d .*log(cosh(s)) .*tanh(s).^2+1512 .*b.^3 .*log(cosh(s)).^2 .*tanh(s).^2);

xi3 = @(s) (1/4802).*a.^(-3).*((-21).*b.*(78.*b.^2+49.*d)+9604.*a.*e+6.*b.*( ...
735.*a1.*b+(-276).*b.^2+490.*((-2).*a.*b1+d)).*log(cosh(s))+7.*( ...
234.*b.^3+147.*b.*d+(-1372).*a.*e+27.*b.*log(cosh(s)).*(6.*b.^2+ ...
49.*d+(-8).*b.^2.*log(cosh(s)))).*sech(s).^2+42.*b.*((-35).*a1.*b+ ...
36.*b.^2+(-98).*d).*s.*tanh(s));

sigma0 = -6;
sigma1 = 0;
sigma2 = 3/49/a^2*(70*a1*b-6*b^2-49*d);
sigma3 = 0;
delta0 = 2;
delta2 = -2*(5*a1*b+7*d)/7/a^2;
omega1 = @(xi) -6/7*b/a*tanh(xi);
omega2 = @(xi) (-72*b^2 + 196*d + 90*b^2*sech(xi).^2 - 441*d*sech(xi).^2 + ...
                70*b*a1 + 144*b^2*tanh(xi).^2)/(196*a^2);
omega3 = @(xi) ((30*a1*b^2)/(49*a^3)-(198*b^3)/(2401*a^3)-(60*b*b1)/(49*a^2)-(33*b*d)/(49*a^3)+(4*e)/a^2)*tanh(xi)+((18*b^3)/(343*a^3)+(3*b*d)/(7*a^3)-(4*e)/a^2)*tanh(xi).^3;
u0  = @(xi) sigma0*sech(xi).^2+delta0;
u2  = @(xi) sigma2*sech(xi).^2+delta2;
v0  = @(xi) 12*sech(xi).^2.*tanh(xi);
v1  = @(xi) (omega1(xi)-sigma1/6).*v0(xi);
v2  = @(xi) (omega2(xi)-sigma1/6*omega1(xi)-sigma2/6).*v0(xi);
v3  = @(xi) (omega3(xi)-sigma2/6*omega1(xi)-sigma1/6*omega2(xi)-sigma3/6).*v0(xi);
switch options.order
				case 1
            xi = @(s) s + xi1(s).*eps;
            u   = @(xi) u0(xi);
            v   = @(xi) v0(xi)+v1(xi).*eps;
        case 2
            xi = @(s) s + xi1(s).*eps + xi2(s).*eps^2;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2;
				otherwise
            xi = @(s) s + xi1(s).*eps + xi2(s).*eps^2 + xi3(s).*eps^3;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2+v3(xi).*eps.^3;
end
if ~HigherOrderTimeReparametrization
    xi = @(s) s;
end
w0  = 1/a*u(xi(eps*t)).*eps.^2;
w1  = 1/a*v(xi(eps*t)).*eps.^3;
ups = q0*w0 + q1*w1 + H0010*beta1 + H0001*beta2 + ...
    1/2*H2000*w0.^2 + H1100*w0.*w1 + 1/2*H0200*w1.^2 +...
    H1010*w0.*beta1 + H1001*w0.*beta2 + H0110*w1.*beta1 + ...
    H0101*w1.*beta2 + 1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
    1/6*H3000*w0.^3 + 1/2*H2100*w0.^2.*w1 + H1101*w0.*w1.*beta2 + ...
    1/2*H2001*w0.^2*beta2 + 1/6*H0003*beta2.^3 + ...
    1/2*H1002*w0*beta2.^2 + 1/2*H0102*w1*beta2.^2;
ups=ups+bt.x; % we shift the orbits with the equilibrium not the saddle!

%% Approximate the sadddle equilibrium
u0inf = delta0;
u2inf = delta2;
uinf  = u0inf+u2inf.*eps.^2;
w0inf = 1/a*uinf.*eps.^2;
homds.x0 = bt.x + q0*w0inf + H0010*beta1 + H0001*beta2 + ...
    1/2*H2000*w0inf.^2  +...
    H1010*w0inf.*beta1 + H1001*w0inf.*beta2 + ...
    1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
    1/6*H3000*w0inf.^3 + ...
    1/2*H2001*w0inf.^2*beta2 + 1/6*H0003*beta2.^3 + ...
    1/2*H1002*w0inf*beta2.^2;

%% Derivative of alpha with respect to epsilon, need to determine the correct
% direction of the tangent vector v
dbeta1depsison = -16/a*eps^3;
dbeta2depsilon = b/a*(2*tau0 + 4*tau2*eps^2)*eps;
dalphadepsilon = K10*dbeta1depsison + K01*dbeta2depsilon + ...
    1/2*K02*2*beta2*dbeta2depsilon + K11*(dbeta1depsison*beta2 + ...
    beta1*dbeta2depsilon) + 1/6*K03*3*beta2^2*dbeta2depsilon;

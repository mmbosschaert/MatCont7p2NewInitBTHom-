function [ups, p, dalphadepsilon] = init_BT_Hom_orbital_with_phase_correction(odefile, bt, ap, options) 

global homds cds

%% Compute normal form coefficients
bt = BT_nmfm_orbital(odefile, bt, ap);
a  = bt.nmfm.a;
b  = bt.nmfm.b;
theta1000 = bt.nmfm.theta1000;
theta0001 = bt.nmfm.theta0001;
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

%% Set amplitude and TTolerance
if ~isempty(options.amplitude)
    amplitude = options.amplitude;
    eps  = sqrt(amplitude*b^2/(6*abs(a)));
else
    eps = options.perturbationparameter;
    amplitude = eps^2/b^2*6*abs(a);
end
if ~isempty(options.TTolerance)
    TTolerance = options.TTolerance;
else
    TTolerance = amplitude*options.amplitudeTToleranceRatio;
end
if options.messages
    fprintf('The initial perturbation parameter epsilon:  %d\n', eps);
    fprintf('The initial amplitude: %g\n', amplitude);
end

%% Approximate parameters
% a) The initical approximation of the parameters
beta1 = -4*a^3/b^4*eps^4;
tau0  = 10/7;
tau2  = 288/2401;

switch options.order
    case 1
        beta2 = a/b*tau0*eps^2;
        alpha = K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2;
    otherwise
        beta2 = a/b*(tau0 + tau2*eps^2)*eps^2;
        alpha = K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2 + K11*beta1.*beta2 + ...
                  1/6*K03*beta2^3;
end
p = bt.par;
p(ap) = p(ap) + alpha;
homds.P0 = p;

%% The initial approximation of cycle
% transformation of time \xi to s
% tau = (2*homds.finemsh - 1) * (homds.T);
% gamma1 =  (1/35).*((-4)+(-59).*(836+15.*4019.^(1/2)).^(-1/3)+(836+15.* ...
%     4019.^(1/2)).^(1/3));
% gamma1 = 0;
% gamma3 = 
% xi1 = @(s) -6.*log(cosh(s))/7;
% xi2 = @(s) (1/98).*((-36).*s+(45+(-21).*gamma1.*(4+7.*gamma1)+72.*log(cosh(s) ...
%     )).*tanh(s));
% xi3 = @(s) (9/2401).*((-105)+(-92).*log(cosh(s))+(91+21.*(3+(-4).*log(cosh( ...
%     s))).*log(cosh(s))).*sech(s).^2+84.*s.*tanh(s));

gamma1 =  (1/35).*((-4)+(-59).*(836+15.*4019.^(1/2)).^(-1/3)+(836+15.* ...
    4019.^(1/2)).^(1/3));
gamma2=0;
gamma3=0;

xi1 = @(s) -6.*log(cosh(s))/7;
xi2 = @(s)(1/98).*((-36).*s+3.*(15+(-7).*gamma1.*( ...
4+7.*gamma1)+24.*log(cosh(s))).*tanh(s));
xi3 = @(s) (1/4802).*((-1890)+98.* ...
gamma1.*(36+63.*gamma1+49.*gamma2.*s)+1656.*log(sech(s))+(-7).*(( ...
-234)+(-7).*gamma1.*((-3)+7.*gamma1).*(9+35.*gamma1)+216.*log( ...
cosh(s)).^2+18.*(9+7.*gamma1.*(4+7.*gamma1)).*log(sech(s))).*sech( ...
s).^2+42.*((-49).*(2+7.*gamma1).*gamma2+36.*s).*tanh(s));

% homoclinic orbit
u0 = @(z) 6*z.^2 - 4; 
u1 = @(z) 0;
u2 = @(z) (18/49).*((-1)+z.^2);
u3 = @(z) (288/3773).*z.*(1+(-1).*z.^2);

v0 = @(z) -12.*z.*(z.^2-1);
v1 = @(z) (72*z.^2.*(-1 + z.^2))/7;
v2 = @(z) (-18*(-1 + z.^2).*(5*z + 9*z.^3))/49;
v3 = @(z) (-72*(-1 + z.^2).*(28 - 678*z.^2 + 231*z.^4))/26411;

% \xi(0)=0 phase condition

%gamma1=0;
%gamma3 = 24/18865;

u1 = @(z) (-12).*gamma1.*z.*((-1)+z.^2);
u2 = @(z) (-6/49).*((-3)+49.*gamma1.^2+98.*gamma2.*z).*((-1)+z.^2);
u3 = @(z) (-12).*gamma3.*z.*((-1)+z.^2);

v1 = @(z) (12/7).*((-1)+z.^2).*((-7).*gamma1+6.* ...
    z.^2+21.*gamma1.*z.^2);
v2 = @(z) (-6/49).*((-1)+z.^2).*(98.*gamma2+15.*z+( ...
    -168).*gamma1.*z+(-245).*gamma1.^2.*z+(-294).*gamma2.*z.^2+27.* ...
    z.^3+336.*gamma1.*z.^3+147.*gamma1.^2.*z.^3);
v3 = @(z) (-6/2401).*((-1)+ ...
    z.^2).*(441.*gamma1+(-4116).*gamma1.^2+(-7203).*gamma1.^3+4802.* ...
    gamma3+(-8232).*gamma2.*z+(-9604).*gamma1.*gamma2.*z+(-648).*z.^2+ ...
    2646.*gamma1.*z.^2+24696.*gamma1.^2.*z.^2+4802.*gamma1.^3.*z.^2+( ...
    -14406).*gamma3.*z.^2+16464.*gamma2.*z.^3+14406.*gamma1.*gamma2.* ...
    z.^3+252.*z.^4+(-6615).*gamma1.*z.^4+(-16464).*gamma1.^2.*z.^4+ ...
    2401.*gamma1.^3.*z.^4);

switch options.order
    case 1
        xi  = @(s) s + xi1(s).*eps;
        u   = @(xi) u0(tanh(xi)) + u1(tanh(xi)).*eps;
        v   = @(xi) v0(tanh(xi)) + v1(tanh(xi)).*eps;
    case 2
        xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2;
        u   = @(xi) u0(tanh(xi))+u1(tanh(xi)).*eps+u2(tanh(xi)).*eps.^2;
        v   = @(xi) v0(tanh(xi))+v1(tanh(xi)).*eps+v2(tanh(xi)).*eps.^2;
    otherwise
        xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2 + xi3(s).*eps^3;
        u   = @(xi) u0(tanh(xi))+u1(tanh(xi)).*eps+u2(tanh(xi)).*eps.^2 ...
            + u3(tanh(xi)).*eps.^3;
        v   = @(xi) v0(tanh(xi))+v1(tanh(xi)).*eps+v2(tanh(xi)).*eps.^2 ...
            + v3(tanh(xi)).*eps.^3;
end
w0  = @(tau)   a/b^2*u(xi(a/b*eps*tau))*eps^2;
w1  = @(tau) a^2/b^3*v(xi(a/b*eps*tau))*eps^3;
%tau = @(t) t./(1 + theta1000.*w0(t) + theta0001.*beta2);

%% tau of t
switch options.order
    case 1
        auxilaryIntegral = @(xi) ...
            (6/7).*eps.*(2.*log(cosh(xi))+(3+(-7).*gamma1).*sech(xi).^2)+ ...
            2.*(xi+(-3).*tanh(xi));

        tOfTau = @(tau) (1+(10/7).*a.*b.^(-1).*eps.^2.*theta0001).*tau ...
            + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
    case 2
        auxilaryIntegral = @(xi) ... 
            (6/7).*eps.*(2.*log(cosh(xi))+(3+(-7).*gamma1).*sech(xi).^2)+ ...
            2.*(xi+(-3).*tanh(xi))+(3/49).*eps.^2.*(12.*xi+(-98).*gamma2.* ...
            sech(xi).^2+((-27)+7.*gamma1.*(4+7.*gamma1)+(15+(-7).*gamma1.*(12+ ...
            7.*gamma1)).*sech(xi).^2).*tanh(xi));

        tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
            + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
    otherwise
        auxilaryIntegral = @(xi) ...
            (6/7).*eps.*(2.*log(cosh(xi))+(3+(-7).*gamma1).*sech(xi).^2)+ ...
            2.*(xi+(-3).*tanh(xi))+(3/49).*eps.^2.*(12.*xi+(-98).*gamma2.* ...
            sech(xi).^2+((-27)+7.*gamma1.*(4+7.*gamma1)+(15+(-7).*gamma1.*(12+ ...
            7.*gamma1)).*sech(xi).^2).*tanh(xi))+(1/2401).*eps.^3.*(144.* ...
            log(cosh(xi))+21.*((-18)+315.*gamma1+343.*gamma1.^3).*sech(xi).^4+ ...
            686.*gamma2.*((-7).*gamma1.*xi+6.*tanh(xi))+(-1).*sech(xi).^2.*(( ...
            -846)+49.*gamma1.*(153+35.*gamma1.*(6+7.*gamma1))+14406.*gamma3+ ...
            2058.*(6+7.*gamma1).*gamma2.*tanh(xi)));

        tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
            + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
end


% Estimate half-return time T
tauofTplus = abs(b/a*1/eps*asech(abs(b)/eps*sqrt(TTolerance/(6*abs(a)))));
homds.T = fzero(@(tau) tauofTplus - tOfTau(tau), tauofTplus);

if options.messages
    fprintf('The initial half-return time T: %g\n', homds.T);
end
t = homds.finemsh;
t = (2*homds.finemsh - 1) * (homds.T);
tauOfT = zeros(size(t));
for i=1:length(t)
    tauOfT(i) = fzero(@(tau) t(i) - tOfTau(tau), t(i));
end

w0 = w0(tauOfT);
w1 = w1(tauOfT);
%w0 = w0(t);
%w1 = w1(t);
ups = q0*w0 + q1*w1 + H0010*beta1 + H0001*beta2 + ...
    1/2*H2000*w0.^2 + H1100*w0.*w1 + 1/2*H0200*w1.^2 +...
    H1010*w0.*beta1 + H1001*w0.*beta2 + H0110*w1.*beta1 + ...
    H0101*w1.*beta2 + 1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
    1/6*H3000*w0.^3 + 1/2*H2100*w0.^2.*w1 + H1101*w0.*w1.*beta2 + ...
    1/2*H2001*w0.^2*beta2 + 1/6*H0003*beta2.^3 + ...
    1/2*H1002*w0*beta2.^2 + 1/2*H0102*w1*beta2.^2;
ups=ups+bt.x; % we shift the orbits with the equilibrium not the saddle!

%% Approximate the sadddle equilibrium
delta0 = 2;
delta2 = 0;
u0inf = delta0;
u2inf = delta2;
uinf  = u0inf+u2inf.*eps.^2;
w0inf = a/b^2*uinf.*eps.^2;
homds.x0 = bt.x + q0*w0inf + H0010*beta1 + H0001*beta2 + ...
    1/2*H2000*w0inf.^2  +...
    H1010*w0inf.*beta1 + H1001*w0inf.*beta2 + ...
    1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
    1/6*H3000*w0inf.^3 + ...
    1/2*H2001*w0inf.^2*beta2 + 1/6*H0003*beta2.^3 + ...
    1/2*H1002*w0inf*beta2.^2;

%% Derivative of alpha with respect to epsilon, need to determine the correct
% direction of the tangent vector v
dbeta1depsison = -16*a^3/b^4*eps^3;
dbeta2depsilon = a/b*(2*tau0 + 4*tau2*eps^2)*eps;
dalphadepsilon = K10*dbeta1depsison + K01*dbeta2depsilon + ...
     1/2*K02*2*beta2*dbeta2depsilon + K11*(dbeta1depsison*beta2 + beta1*dbeta2depsilon) + ...
     1/6*K03*3*beta2^2*dbeta2depsilon;

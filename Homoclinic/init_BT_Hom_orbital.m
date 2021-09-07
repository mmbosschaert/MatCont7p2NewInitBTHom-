function [ups, p, dalphadepsilon] = init_BT_Hom_orbital(odefile, bt, ap, options) 

global homds cds

%% Compute normal form coefficients
bt = BT_nmfm_orbital(odefile, bt, ap, options);
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
if ~isempty(options.amplitude) && options.amplitude ~= 0
    amplitude = options.amplitude;
    eps = sqrt(amplitude*b^2/(6*abs(a)));
else
    eps = options.perturbationparameter;
    amplitude = eps^2/b^2*6*abs(a);
end
if ~isempty(options.TTolerance) && options.TTolerance ~= 0
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
xi1 = @(s) -6.*log(cosh(s))/7;
xi2 = @(s) (-9/98).*(4.*s+(-5).*tanh(s)+(-8).*log( cosh(s)).*tanh(s));
xi3 = @(s) (9/2401).*((-29).*log(cosh(s))+(-84).*log( ...
cosh(s)).^2+(-14).*sech(s).^2+84.*s.*tanh(s)+(-105).*tanh(s).^2+( ...
-63).*log(cosh(s)).*tanh(s).^2+84.*log(cosh(s)).^2.*tanh(s).^2);

xi3 = @(s) (9/2401).*((-91)+(-92).*log(cosh(s))+(91+21.*(3+(-4).*log(cosh(s)) ...
    ).*log(cosh(s))).*sech(s).^2+84.*s.*tanh(s));

% homoclinic orbit
u0 = @(xi) 6*tanh(xi).^2 - 4; 
u2 = @(xi) -18/49*sech(xi).^2; 
v0 = @(xi) 12.*sech(xi).^2.*tanh(xi);
v1 = @(xi) -72/7*tanh(xi).*sech(xi).^2.*tanh(xi);
v2 = @(xi) (90/49 + 162/49*tanh(xi).^2).*sech(xi).^2.*tanh(xi);
v3 = @(xi) -(3888/2401*tanh(xi) - 216/343*tanh(xi).^3).*sech(xi).^2.*tanh(xi);
switch options.order
				case 1
            xi  = @(s) s + xi1(s).*eps;
            u   = @(xi) u0(xi);
            v   = @(xi) v0(xi)+v1(xi).*eps;
        case 2
            xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2;
				otherwise
            xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2 + xi3(s).*eps^3;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2+v3(xi).*eps.^3;
end
w0  = @(tau)   a/b^2*u(xi(a/b*eps*tau))*eps^2;
w1  = @(tau) a^2/b^3*v(xi(a/b*eps*tau))*eps^3;

%% tau of t
switch options.order
		case 1
        auxilaryIntegral = @(xi) ...
            2.*xi-6*tanh(xi)+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2);

        tOfTau = @(tau) (1+(10/7).*a.*b.^(-1).*eps.^2.*theta0001).*tau ...
                    + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
    case 2
        auxilaryIntegral = @(xi) ... 
            2.*xi+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2)+( ...
            -6).*tanh(xi)+eps.^2.*((36/49).*xi+(-81/49).*tanh(xi)+(45/49).* ...
            sech(xi).^2.*tanh(xi));

        tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
                    + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
		otherwise
        auxilaryIntegral = @(xi) ...
            2.*xi+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2)+ ...
            eps.^3.*((-468/2401)+(144/2401).*log(cosh(xi))+(846/2401).*sech( ...
            xi).^2+(-54/343).*sech(xi).^4)+(-6).*tanh(xi)+eps.^2.*((36/49).* ...
            xi+(-81/49).*tanh(xi)+(45/49).*sech(xi).^2.*tanh(xi));

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

% filename = 'tauOfT.csv';
% writematrix([t; tauOfT]', filename, 'Delimiter',  ' ');
% tau = linspace(min(tauOfT),max(tauOfT),60);
% filename = 'tOfTau.csv';
% writematrix([tau; tOfTau(tau)]', filename, 'Delimiter',  ' ');
% figure; hold on
% plot(t, tauOfT)
% plot(tau, tOfTau(tau))

w0 = w0(tauOfT);
w1 = w1(tauOfT);
switch options.order
    case 1
        ups = q0*w0 + q1*w1 + H0010*beta1 + H0001*beta2 + ...
            1/2*H2000*w0.^2 + H1100*w0.*w1 + ...
            H1001*w0.*beta2 + H0101*w1.*beta2 + 1/2*H0002*beta2.^2;
    case 2
        ups = q0*w0 + q1*w1 + H0010*beta1 + H0001*beta2 + ...
            1/2*H2000*w0.^2 + H1100*w0.*w1 + 1/2*H0200*w1.^2 +...
            H1010*w0.*beta1 + H1001*w0.*beta2 +  ...
            H0101*w1.*beta2 + 1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
            1/6*H3000*w0.^3 + ...
            1/2*H2001*w0.^2*beta2 + 1/6*H0003*beta2.^3 + ...
            1/2*H1002*w0*beta2.^2;
    otherwise
        ups = q0*w0 + q1*w1 + H0010*beta1 + H0001*beta2 + ...
            1/2*H2000*w0.^2 + H1100*w0.*w1 + 1/2*H0200*w1.^2 +...
            H1010*w0.*beta1 + H1001*w0.*beta2 + H0110*w1.*beta1 + ...
            H0101*w1.*beta2 + 1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
            1/6*H3000*w0.^3 + 1/2*H2100*w0.^2.*w1 + H1101*w0.*w1.*beta2 + ...
            1/2*H2001*w0.^2*beta2 + 1/6*H0003*beta2.^3 + ...
            1/2*H1002*w0*beta2.^2 + 1/2*H0102*w1*beta2.^2;
end
ups=ups+bt.x; % we shift the orbits with the equilibrium not the saddle!

%% Approximate the sadddle equilibrium
delta0 = 2;
delta2 = 0;
u0inf = delta0;
u2inf = delta2;
uinf  = u0inf+u2inf.*eps.^2;
w0inf = a/b^2*uinf.*eps.^2;
switch options.order
    case 1
        homds.x0 = bt.x + q0*w0inf + H0010*beta1 + H0001*beta2 + ...
            1/2*H2000*w0inf.^2 + H1001*w0inf.*beta2 + 1/2*H0002*beta2.^2;
    otherwise
        homds.x0 = bt.x + q0*w0inf + H0010*beta1 + H0001*beta2 + ...
            1/2*H2000*w0inf.^2  +...
            H1010*w0inf.*beta1 + H1001*w0inf.*beta2 + ...
            1/2*H0002*beta2.^2 + H0011*beta1.*beta2 + ...
            1/6*H3000*w0inf.^3 + ...
            1/2*H2001*w0inf.^2*beta2 + 1/6*H0003*beta2.^3 + ...
            1/2*H1002*w0inf*beta2.^2;
end

%% Derivative of alpha with respect to epsilon, need to determine the correct
% direction of the tangent vector v
dbeta1depsilon = -16*a^3/b^4*eps^3;
dbeta2depsilon = a/b*(2*tau0 + 4*tau2*eps^2)*eps;
dalphadepsilon = K10*dbeta1depsilon + K01*dbeta2depsilon + ...
     K02*beta2*dbeta2depsilon + K11*(dbeta1depsilon*beta2 + beta1*dbeta2depsilon) + ...
     1/6*K03*3*beta2^2*dbeta2depsilon;

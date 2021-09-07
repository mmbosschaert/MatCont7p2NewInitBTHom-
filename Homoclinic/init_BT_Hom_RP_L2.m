function [ups, p, dalphadepsilon] = init_BT_Hom_RP_L2(odefile, bt, ap, options)

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

%% Approximation of homoclic solution
% transformation of time \xi to s
u0 = @(s) 2*(1-3*sech(s).^2);
v0 = @(s) 2*(6*sech(s).^2.*tanh(s));

u1 = @(s) (36/245).*(59+(-70).*log(2.*cosh(s))).*sech(s).^2.*tanh(s);
v1 = @(s) (36/245).*(153+(-140).*log(2.*cosh(s))+cosh(2.*s).*((-94)+70.*log( ...
2.*cosh(s)))).*sech(s).^4;

u2 = @(s) (36/60025).*sech(s).^2.*(3.*(6289+70.*log(2.*cosh(s)).*((-247)+ ...
105.*log(2.*cosh(s)))).*sech(s).^2+(-2).*(7129+210.*log(2.*cosh(s) ...
).*((-94)+35.*log(2.*cosh(s)))+3675.*s.*tanh(s)));
v2 = @(s) (18/60025).*sech(s).^5.*((-22050).*s.*cosh(s)+7350.*s.*cosh(3.*s)+ ...
2.*((-97015)+197400.*log(2.*cosh(s))+(-73500).*log(2.*cosh(s)).^2+ ...
cosh(2.*s).*(30323+(-54180).*log(2.*cosh(s))+14700.*log(2.*cosh(s) ...
).^2)).*sinh(s));

u3 = @(s) (2/1191196125).*sech(s).^2.*(8748.*((-70).*(210.*s.*((-47)+35.* ...
log(2.*cosh(s)))+30673.*log(cosh(s)).*tanh(s))+sech(s).^2.*(3675.* ...
s.*((-247)+210.*log(2.*cosh(s)))+((-966242)+4456830.*log(2.*cosh( ...
s))+7350.*((-470)+129.*cosh(2.*s)).*log(2.*cosh(s)).^2+(-171500).* ...
((-5)+cosh(2.*s)).*log(2.*cosh(s)).^3).*tanh(s)))+(-1).*tanh(s).*( ...
884895199+231525000.*pi.^2+18782918280.*log(2)+(-5501034000).* ...
zeta(3)));

v3 = @(s) (1/2382392250).*sech(s).^6.*(cosh(4.*s).*((-5484567341)+ ...
231525000.*pi .^2+18782918280.*log(2)+18782918280.*log(cosh(s))+ ...
14338409400.*log(2.*cosh(s))+(-21089678400).*log(2.*cosh(s)).^2+ ...
3000564000.*log(2.*cosh(s)).^3+(-5501034000).*zeta(3))+(-3).*( ...
59829237043+231525000.*pi .^2+18782918280.*log(2)+18782918280.* ...
log(cosh(s))+(-204026717160).*log(2.*cosh(s))+153971798400.*log( ...
2.*cosh(s)).^2+(-33006204000).*log(2.*cosh(s)).^3+428652000.*s.*(( ...
-47)+35.*log(2.*cosh(s))).*sinh(2.*s)+(-21432600).*s.*((-129)+70.* ...
log(2.*cosh(s))).*sinh(4.*s)+(-5501034000).*zeta(3))+(-10).*cosh( ...
2.*s).*((-14003481721)+46305000.*pi .^2+3756583656.*log(2)+ ...
3756583656.*log(cosh(s))+46146347352.*log(2.*cosh(s))+( ...
-39530287440).*log(2.*cosh(s)).^2+7801466400.*log(2.*cosh(s)).^3+( ...
-1100206800).*zeta(3)));

switch options.order
  case 0
    u  = @(s) u0(s);
    v  = @(s) v0(s);
  case 1
    u = @(s) u0(s) + eps.*u1(s);
    v = @(s) v0(s) + eps.*v1(s);
  case 2
    u = @(s) u0(s) + eps.*u1(s) + eps^2.*u2(s);
    v = @(s) v0(s) + eps.*v1(s) + eps^2.*v2(s);
  case 3
    u = @(s) u0(s) + eps.*u1(s) + eps^2.*u2(s) + eps^3.*u3(s);
    v = @(s) v0(s) + eps.*v1(s) + eps^2.*v2(s) + eps^3.*v3(s);
end

w0  = @(tau)   a/b^2*u(a/b*eps*tau)*eps^2;
w1  = @(tau) a^2/b^3*v(a/b*eps*tau)*eps^3;

%% tau of t
switch options.order
		case 1
        auxilaryIntegral = @(s) ...
            2.*s+(-36/245).*eps.*((-12)+35.*log(2)+(12+(-35).*log(2.*cosh(s))) ...
            .*sech(s).^2)+(-6).*tanh(s);

        tOfTau = @(tau) (1+(10/7).*a.*b.^(-1).*eps.^2.*theta0001).*tau ...
                    + theta1000*1/b*eps*auxilaryIntegral(a/b*tau*eps);
    case 2
        auxilaryIntegral = @(s) ... 
            2.*s+(-36/245).*eps.*((-12)+35.*log(2)+(12+(-35).*log(2.*cosh(s))) ...
            .*sech(s).^2)+(-6).*tanh(s)+(9/60025).*eps.^2.*sech(s).^2.*( ...
            14700.*s+(-1225).*sech(s).*sinh(3.*s)+(7411+(-49560).*log(2.*cosh( ...
            s))+29400.*log(2.*cosh(s)).^2).*tanh(s));

        tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
                    + theta1000*1/b*eps*auxilaryIntegral(a/b*tau*eps);
		otherwise
        zeta3 = zeta(3); % calling zeta(3) each time takes a very long time
        auxilaryIntegral = @(s) ...
            2.*s+(-36/245).*eps.*((-12)+35.*log(2)+(12+(-35).*log(2.*cosh(s))) ...
            .*sech(s).^2)+(-6).*tanh(s)+(9/60025).*eps.^2.*sech(s).^2.*( ...
            14700.*s+(-1225).*sech(s).*sinh(3.*s)+(7411+840.*log(2.*cosh(s)).* ...
            ((-59)+35.*log(2.*cosh(s)))).*tanh(s))+(1/1191196125).*eps.^3.*(( ...
            -4798781017)+(-231525000).*pi .^2+612360.*log(2).*(4259+35.*log(2) ...
            .*((-177)+70.*log(2)))+(-4374).*((-305897)+420.*log(2.*cosh(s)).*( ...
            6289+35.*log(2.*cosh(s)).*((-247)+70.*log(2.*cosh(s))))).*sech(s) ...
            .^4+sech(s).^2.*(3460787539+231525000.*pi .^2+2449440.*log(2.* ...
            cosh(s)).*(3652+35.*log(2.*cosh(s)).*((-141)+35.*log(2.*cosh(s)))) ...
            +64297800.*s.*((-59)+70.*log(2.*cosh(s))).*tanh(s)+(-5501034000).* ...
            zeta3)+5501034000.*zeta3);

        tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
                    + theta1000*1/b*eps*auxilaryIntegral(a/b*tau*eps);
end

% Approximate half-return time T
tauofTplus = abs(b/a*1/eps*asech(abs(b)/eps*sqrt(TTolerance/(6*abs(a)))));
homds.T = fzero(@(tau) tauofTplus - tOfTau(tau), tauofTplus);

if options.messages
    fprintf('The initial perturbation parameter:  %d\n', eps);
    fprintf('The initial amplitude: %g\n', amplitude);
    fprintf('T: %g\n', homds.T);
end
t = homds.finemsh;
t = (2*homds.finemsh - 1) * (homds.T);
tauOfT = zeros(size(t));
dtdtau = @(tau) -(1+theta0001*beta2 + theta1000*w0(tau));
for i=1:length(t)
    tauOfT(i) = fzero(@(tau) t(i) - tOfTau(tau), t(i));
end

% tau = linspace(min(tauOfT),max(tauOfT),60);
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
dbeta1depsison = -16*a^3/b^4*eps^3;
dbeta2depsilon = a/b*(2*tau0 + 4*tau2*eps^2)*eps;
dalphadepsilon = K10*dbeta1depsison + K01*dbeta2depsilon + ...
     1/2*K02*2*beta2*dbeta2depsilon + K11*(dbeta1depsison*beta2 + beta1*dbeta2depsilon) + ...
     1/6*K03*3*beta2^2*dbeta2depsilon;

matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'HodgkinHuxley';

coordsnames = {'V', 'm', 'n', 'h'};
parnames={'VbarK', 'I'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

gbarNa = 120;
gbarK  = 36;
gbarL  = 0.3;
VbarNa = -115;
VbarL  = 10.599;
T = 6.3;

Psi = @(x) x/(exp(x)-1);
alpha_m = @(V) Psi( (V+25)/10 );
alpha_n = @(V) 0.1*Psi( (V+10)/10);
alpha_h = @(V) 0.07*exp(V/20);

beta_m = @(V) 4*exp(V/18);
beta_n = @(V) 0.125*exp(V/80);
beta_h = @(V) 1/(1+exp((V+30)/10));

G = @(V, m, n, h) gbarNa*m^3*h*(V-VbarNa) + gbarK*n^4*(V-VbarK) + gbarL*(V-VbarL);
Phi = @(T) 3^(T-6.3)/10;
dV_dt = -G(V, m, n, h)+I;
dm_dt = Phi(T)*((1-m)*alpha_m(V)-m*beta_m(V));
dn_dt = Phi(T)*((1-n)*alpha_n(V)-n*beta_n(V));
dh_dt = Phi(T)*((1-h)*alpha_h(V)-h*beta_h(V));
system = [dV_dt; dm_dt; dn_dt; dh_dt];

suc = generate_directional_derivatives(...
  system,...   % n x 1 array of derivative symbolic expressions
  coords,... % 1 x n array of symbols for states
  par,...      % 1 x np array of symbols used for parameters
  system_name,... % argument specifying the system name
  [matcontpath, 'Systems/']... % directory to save to file
);

order = 3;
suc = generate_multilinear_forms(system_name, system, coords, par, order, ...
        [matcontpath, 'Systems/']);

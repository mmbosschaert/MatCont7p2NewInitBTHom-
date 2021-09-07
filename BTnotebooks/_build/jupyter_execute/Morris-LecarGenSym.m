matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'Morris_Lecar';

coordsnames = {'V', 'w'};
parnames = {'Iapp', 'v3'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

C=20;
VL=-60;
VCa=120;
VK=-84;
gL=2;
gCa=4.4;
gK=8;
v1=-1.2;
v2=18;
v4=30;
phi=1/25;

minf = 0.5*(1+tanh((V-v1)/v2));
winf = 0.5*(1+tanh((V-v3)/v4));
Iion = gCa*minf*(V-VCa) + gK*w*(V-VK) + gL*(V-VL);
tau = sech((V-v3)/2/v4);
dV_dt = (Iapp - Iion)/C;
dw_dt = phi/tau*(winf-w);
system = [dV_dt; dw_dt];

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

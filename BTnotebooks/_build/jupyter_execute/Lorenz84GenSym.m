matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, 'Utilities'])
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'extendedLorenz84';

coordsnames = {'X', 'Y', 'Z', 'U'};
parnames = {'F', 'S'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

alpha = 0.25;
beta = 1;
G = 0.25;
delta = 1.04;
xi = 0.987;

dX_dt = -Y^2-Z^2-alpha*X+alpha*F-xi*U^2;
dY_dt = X*Y-beta*X*Z-Y+G;
dZ_dt = beta*X*Y+X*Z-Z;
dU_dt = -delta*U+xi*U*X+S;
system = [dX_dt; dY_dt; dZ_dt; dU_dt];

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

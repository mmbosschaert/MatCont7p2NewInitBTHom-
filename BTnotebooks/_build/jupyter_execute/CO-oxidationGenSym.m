matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'CO_oxidation';

coordsnames = {'x', 'y', 's'};
parnames =  {'k1', 'km1', 'k3', 'k2', 'km2', 'k4', 'lambda'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

z = 1-x-y-s;
dx_dt = 2*k1*z^2 - 2*km1*x^2 - k3*x*y;
dy_dt = k2*z - km2*y - k3*x*y;
ds_dt = k4*(z - lambda*s);
system = [dx_dt; dy_dt; ds_dt];

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

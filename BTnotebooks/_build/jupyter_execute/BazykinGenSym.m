matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'Bazykin';

coordsnames = {'x1', 'x2'};
parnames = {'alpha', 'delta'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

gamma = 1;
epsilon = 0.01;

dx1_dt = x1 - x1*x2/(1+alpha*x1) - epsilon*x1^2;
dx2_dt = -gamma*x2 + x1*x2/(1+alpha*x1) - delta*x2^2;
system = [dx1_dt; dx2_dt];

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

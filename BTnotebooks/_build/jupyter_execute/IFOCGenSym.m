matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end

system_name = 'IFOC';

coordsnames = {'x1', 'x2', 'x3', 'x4'};
parnames={'k', 'Tm'};

syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates

c1 = 4.4868;
c2 = 0.3567;
c3 = 0;
c4 = 9.743;
c5 = 1.911;
u20 = 11.3;
kp = 4.5;
ki = 500;
wref = 0;

dx1_dt = -c1*x1 + c2*x4 - k*c1/u20*x2*x4;
dx2_dt = -c1*x2 + c2*u20 + k*c1/u20*x1*x4;
dx3_dt = -c3*x3 - c4*c5*(x2*x4 - u20*x1) + (c4*Tm + c3*wref);
dx4_dt = -(ki-kp*c3)*x3 - kp*c4*c5*(x2*x4 - u20*x1) + kp*(c4*Tm + c3*wref);
system = [dx1_dt; dx2_dt; dx3_dt; dx4_dt];

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

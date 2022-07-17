clear all
matcontpath = '../';
addpath(matcontpath)
addpath([matcontpath, 'Equilibrium'])
addpath([matcontpath, 'Systems'])
addpath([matcontpath, 'Hopf'])
addpath([matcontpath, 'Homoclinic'])
addpath([matcontpath, 'LimitPoint'])
addpath([matcontpath, 'LimitCycle'])
addpath([matcontpath, 'Continuer'])
addpath([matcontpath, 'MultilinearForms'])
addpath([matcontpath, 'Utilities'])
set(groot, 'defaultTextInterpreter', 'LaTeX');
set(0,'defaultAxesFontSize',15)

odefile=@SIRmodel;

A     =  20;
mu0   =  10;
d     =  1/10;
nu    =  1;
beta  =  11.5;

mu  =  @(b,mu1,I) mu0 + (mu1-mu0)*b/(I+b);
S   =  @(b,I,mu1) (A-(d+nu+mu(b, mu1, I))*I)/d; 
R   =  @(b,I,mu1)   mu(b, mu1, I)*I/d;

R0 = @(mu1) beta/(d+nu+mu1);
Ascr = (d+nu+mu0)*(beta-nu);
Bscr = @(b,mu1) (d+nu+mu0-beta)*A+(beta-nu)*(d+nu+mu1)*b;
Cscr = @(b,mu1) (d+nu+mu1)*A*b*(1-R0(mu1));

delta0 = d+nu+mu0;
delta1 = @(mu1) d+nu+mu1;
Delta0 = @(b,mu1) (beta-nu)^2*delta1(mu1)^2*b^2-2*A*(beta-nu)*(beta*(mu1-mu0)+delta0*(delta1(mu1)-beta))*b+A^2*(beta-delta0)^2;
I1 = @(b,mu1)(-Bscr(b,mu1)-sqrt(Delta0(b,mu1)))/(2*Ascr);
I2 = @(b,mu1)(-Bscr(b,mu1)+sqrt(Delta0(b,mu1)))/(2*Ascr);

b=0.022;
mu1=10.45;
p(1) = mu1;
p(2) = b;
x = [S(b,I2(b,mu1),mu1); I2(b,mu1); R(b,I2(b,mu1),mu1)];

parnames = {'mu1','b'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});
ap = [ind.mu1, ind.b]; % continuation parameters

[x1_pred, v1_pred] = init_EP_EP(odefile, x, p, ind.b);
opt = contset;
opt.MaxNumPoints  = 1000;
opt.Singularities = 1;
opt.TestTolerance = 1e-15;
opt.MaxTestIters = 10;
opt.MinStepsize = 0.001;
[eqbr_x, ~, eqbr_bif_data] = cont(@equilibrium, x1_pred, v1_pred, opt);

%plot --width 1024 --height 800
hopfPointInfo   = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Hopf')==1);
limitPointsInfo = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Limit point')==1);
HopfPoint1 = eqbr_x(:,hopfPointInfo.index);
limitpoint1 = eqbr_x(:,limitPointsInfo.index);
plot(eqbr_x(4,:), eqbr_x(1,:)); hold on
plot(HopfPoint1(4), HopfPoint1(1), '.r' ,'MarkerSize', 20)
plot(limitpoint1(4), limitpoint1(1),  '.k' ,'MarkerSize', 20)
xlabel('$b$')
ylabel('$S$')
legend({'Equilibrium curve', 'Hopf point', 'Limit point'}, 'Location', 'NorthWest')
title('Equilibrium curve in $(b,S)$-space', 'fontsize', 15)

hopf1.par = [mu1, b];
hopf1.par(ind.b) = eqbr_x(4, hopfPointInfo(1).index);
hopf1.x = eqbr_x(1:3, hopfPointInfo(1).index);
hopf1.x = x;
[hopf_x, hopf_v] = init_H_H(odefile, hopf1.x, hopf1.par', ap);
opt.MaxNumPoints = 500;
opt.MaxStepsize = .01;
opt.Singularities = 1;
[hopf_br, ~, hopf_br_bif] = cont(@hopf, hopf_x, hopf_v, opt);

bt_points_info = hopf_br_bif(strcmp({hopf_br_bif.label}, 'BT')==1);
bt.x = hopf_br(1:3, bt_points_info.index);
bt.par = hopf_br(4:5, bt_points_info.index);
BToptions = BT_Hom_set_options();
[hom_x, hom_v] = init_BT_Hom(odefile, bt,  ap, BToptions);
opt = contset;
opt.Singularities = 0;
homoclinic_br = cont(@homoclinic, hom_x, hom_v, opt);

global homds
hold on
% plot computed homoclinic parameter curve
plot(homoclinic_br(homds.PeriodIdx+1,:), ...
     homoclinic_br(homds.PeriodIdx+2,:));
% Bogdanov-Takens parameter-dependent smooth orbital normal form coefficients
bt = BT_nmfm_orbital(odefile, bt, ap, BToptions);
a   = bt.nmfm.a;
b   = bt.nmfm.b;
K10 = bt.nmfm.K10;
K01 = bt.nmfm.K01;
K02 = bt.nmfm.K02;
K11 = bt.nmfm.K11;
K03 = bt.nmfm.K03;
% construct predictor as in the paper
eps = linspace(0, 0.3);
beta1 = -4*a^3/b^4*eps.^4;
tau0  = 10/7;
tau2  = 288/2401;
beta2 = a/b*(tau0 + tau2*eps.^2).*eps.^2;
alpha = K10.*beta1 + K01.*beta2 + 1/2*K02.*beta2.^2 ...
    + K11.*beta1.*beta2 + 1/6*K03.*beta2.^3;
alpha = bt.par(ap) + alpha;
% plot currect predictor
plot(alpha(1,:), alpha(2,:), '.')
% plot Bogdanov-Takens point
plot(bt.par(ind.mu1), bt.par(ind.b), '.k', 'MarkerSize', 20)
plot(hopf_br(4,:), hopf_br(5,:))
% set axis labels and legend
xlabel('$\mu_1$')
ylabel('$b$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point', 'Hopf/Neutral saddle curve'}, 'Location', 'SouthWest')
title('Comparision between computed and predicted parameter curve.')

hold on
plot3(homoclinic_br(homds.coords(1:homds.nphase:end), 1:10:end), ...
      homoclinic_br(homds.coords(2:homds.nphase:end), 1:10:end), ...
      homoclinic_br(homds.coords(3:homds.nphase:end), 1:10:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$S$')
ylabel('$I$')
zlabel('$R$')
plot3(bt.x(1), bt.x(2), bt.x(3), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(x_1,x_2,x_3)$ phase space')
grid on
view([73, 33]);

options = BT_Hom_set_options();
options.messages = false;
options.correct = false;
options.TTolerance = 1.0e-05;

amplitudes = linspace(1.0e-03, 1,10);
XPredicted = zeros(494,length(amplitudes));
XCorrected = zeros(494,length(amplitudes));
for j=1:length(amplitudes)
  options.amplitude = amplitudes(j);
  [x_pred, v0] = init_BT_Hom(odefile, bt, ap, options);
  XPredicted(:,j) = x_pred;
  try
    XCorrected(:,j) = newtcorr(x_pred, v0);
  catch
    warning('Didn''t convergence to homoclinic solution')
  end
end

hold on
cm = lines;
plot(XPredicted(homds.coords(1:homds.nphase:end),1:10), ...
     XPredicted(homds.coords(2:homds.nphase:end),1:10), ...
      'color', cm(1,:))
plot(XCorrected(homds.coords(1:homds.nphase:end),1:10), ...
     XCorrected(homds.coords(2:homds.nphase:end),1:10), ...
      '--', 'color', cm(2,:))
plot(bt.x(1), bt.x(2), '.', 'MarkerSize', 16)
xlabel('$S$')
ylabel('$I$')
zlabel('$R$')
grid on
view([73, 33]);

BToptions = BT_Hom_set_options();
BToptions.TTolerance = 1e-05;
BToptions.messages = false;
BToptions.correct = false;

amplitudes = logspace(-3, 1, 20);
methodList = {'orbital', 'LP', 'RegularPerturbation', ...
    'RegularPerturbationL2', 'LPHypernormalForm'};
relativeErrors = {};
for i=1:length(methodList)
    BToptions.method = methodList{i};
    relativeErrors{i} = zeros(size(amplitudes));
    for j=1:length(amplitudes)
    BToptions.amplitude = amplitudes(j);
    [x_pred, v0] = init_BT_Hom(odefile, bt, ap, BToptions);
    try
        x_corrected = newtcorr(x_pred, v0);
        relativeErrors{i}(j) = norm(x_corrected-x_pred)/norm(x_corrected);
    catch
        warning('Did not converge.')
        continue
    end
  end
end

cm = lines();
loglog(amplitudes, relativeErrors{1}(:), 'd', ...
       amplitudes, relativeErrors{2}(:), '--', ...
       amplitudes, relativeErrors{3}(:), '*', ...
       amplitudes, relativeErrors{4}(:), 's', ...
       amplitudes, relativeErrors{5}(:), '+')
legend(methodList, 'Location', 'NorthWest')
title('SIR model')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
ax.ColorOrder = [cm(1,:); [0.8 0.8 0.8]; cm(2,:); cm(4,:); cm(5,:)];

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

odefile=@IFOC;

c1 = 4.4868;
c2 = 0.3567;
c5 = 1.911;
u20 = 11.3;
k = 1+(1/3)*(216+(-24)*57^(1/2))^(1/3)+2*3^(-2/3)*(9+57^(1/2))^(1/3);
x1m = (1/4)*2^(-1/2)*c1^(-1)*c2*(1+k)^(-1)*((-1)+k^2+(-1)*(9+( ...
      -10)*k^2+k^4)^(1/2))*(k^(-2)*((-3)+k^2+(9+(-10)*k^2+ ...
      k^4)^(1/2))*u20^2)^(1/2);
x2m = (1/4)*c1^(-1)*c2*k^(-1)*(1+k)^(-1)*(3+k*(4+k)+(-1)*(9+( ...
      -10)*k^2+k^4)^(1/2))*u20;
x3 = 0;
x4p = (-1)*2^(-1/2)*(k^(-2)*((-3)+k^2+(9+(-10)*k^2+k^4)^(1/2))*u20^2)^(1/2);
Tmm = (-1/4)*2^(-1/2)*c1^(-1)*c2*c5*k^(-1)*(3+k^2+(-1)*(9+( ...
      -10)*k^2+k^4)^(1/2))*u20*(k^(-2)*((-3)+k^2+(9+(-10)* ...
      k^2+k^4)^(1/2))*u20^2)^(1/2);
bt1.x = [x1m; x2m; x3; x4p];
bt1.par = [k; Tmm];

parnames = {'k', 'Tm'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});

ap = [ind.k, ind.Tm];
BToptions = BT_Hom_set_options();
[hom_x, hom_v] = init_BT_Hom(odefile, bt1,  ap, BToptions);
opt = contset;
opt.Singularities = 0;
opt.MaxNumPoints = 900;
homoclinic_br1 = cont(@homoclinic, hom_x, hom_v, opt);

%plot --width 1024 --height 800
hold on
global homds
% plot computed parameter curve
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:));
% Bogdanov-Takens parameter-dependent normal form coefficients
bt1 = BT_nmfm_orbital(odefile, bt1, ap, BToptions);
a   = bt1.nmfm.a;
b   = bt1.nmfm.b;
K10 = bt1.nmfm.K10;
K01 = bt1.nmfm.K01;
K02 = bt1.nmfm.K02;
K11 = bt1.nmfm.K11;
K03 = bt1.nmfm.K03;
% construct predictor as in the paper
eps = linspace(0, 0.15, 30);
beta1 = -4*a^3/b^4*eps.^4;
tau0  = 10/7;
tau2  = 288/2401;
beta2 = a/b*(tau0 + tau2*eps.^2).*eps.^2;
alpha = K10.*beta1 + K01.*beta2 + 1/2*K02.*beta2.^2 ...
    + K11.*beta1.*beta2 + 1/6*K03.*beta2.^3;
alpha = bt1.par(ap) + alpha;
% plot currect predictor
plot(alpha(1,:), alpha(2,:), '.')
% plot Bogdanov-Takens point
plot(bt1.par(ind.k), bt1.par(ind.Tm), '.k', 'MarkerSize', 20)
% set axis labels and legend
xlabel('$I_app$')
ylabel('$h$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point'}, 'Location', 'SouthWest')
title('Comparision between computed and predicted parameter curve.')

[hom_x, hom_v] = init_BT_Hom(odefile, bt1,  ap, BToptions);
opt.Singularities = 0;
opt.MaxNumPoints = 120;
opt.Backward = 1;
homoclinic_br1_rev = cont(@homoclinic, hom_x, hom_v, opt);

%plot --width 1024 --height 800
hold on
cm = lines();
homColor  = cm(1,:);
% plot computed parameter curve
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:), ...
     'Color', homColor);
plot(homoclinic_br1_rev(homds.PeriodIdx+1,:), ...
     homoclinic_br1_rev(homds.PeriodIdx+2,:), ...
     'Color', homColor);
% plot currect predictor
plot(alpha(1,:), alpha(2,:), '.', 'Color', cm(2,:))
% plot Bogdanov-Takens point
plot(bt1.par(ind.k), bt1.par(ind.Tm), '.k', 'MarkerSize', 20)
% set axis labels and legend
xlabel('$I_{app}$')
ylabel('$h$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point'}, 'Location', 'SouthWest')
title('Comparision between computed and predicted parameter curve.')

hold on
plot3(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:30:end), ...
      homoclinic_br1(homds.coords(2:homds.nphase:end), 1:30:end), ...
      homoclinic_br1(homds.coords(3:homds.nphase:end), 1:30:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$x_3$')
plot3(bt1.x(1), bt1.x(2), bt1.x(3), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(x_1,x_2,x_3)$ phase space')
grid on
view([144, 31]);

options = BT_Hom_set_options();
options.messages = false;
options.correct = false;
options.TTolerance = 1.0e-05;

amplitudes = linspace(1.0e-03, 1, 10);
XPredicted = zeros(660,length(amplitudes));
XCorrected = zeros(660,length(amplitudes));
for j=1:length(amplitudes)
  options.amplitude = amplitudes(j);
  [x_pred, v0] = init_BT_Hom(odefile, bt1, ap, options);
  XPredicted(:,j) = x_pred;
  try
    XCorrected(:,j) = newtcorr(x_pred, v0);
  catch
    warning('Didn''t convergence to homoclinic solution')
  end
end

hold on
cm = lines;
plot3(XPredicted(homds.coords(1:homds.nphase:end),1:10), ...
      XPredicted(homds.coords(2:homds.nphase:end),1:10), ...
      XPredicted(homds.coords(3:homds.nphase:end),1:10), ...
      'color', cm(1,:), 'HandleVisibility', 'Off')
plot3(XCorrected(homds.coords(1:homds.nphase:end),1:10), ...
      XCorrected(homds.coords(2:homds.nphase:end),1:10), ...
      XCorrected(homds.coords(3:homds.nphase:end),1:10), ...
      '--', 'color', cm(2,:), 'HandleVisibility', 'Off')
plot3(bt1.x(1), bt1.x(2), bt1.x(3), '.', 'MarkerSize', 16)
legend('Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(x_1,x_2,x_3)$ phase space')
xlabel('$x_1$')
xlabel('$x_2$')
xlabel('$x_3$')
grid on
view([160, 45]);

BToptions = BT_Hom_set_options();
BToptions.TTolerance = 1e-05;
BToptions.messages = false;
BToptions.correct = false;

amplitudes = logspace(-4, 1, 20);
methodList = {'orbital', 'LP', 'RegularPerturbation', ...
    'RegularPerturbationL2', 'LPHypernormalForm'};
relativeErrors = {};
for i=1:length(methodList)
    BToptions.method = methodList{i};
    relativeErrors{i} = zeros(size(amplitudes));
    for j=1:length(amplitudes)
    BToptions.amplitude = amplitudes(j);
    [x_pred, v0] = init_BT_Hom(odefile, bt1, ap, BToptions);
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
title('IFOC model')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
ax.ColorOrder = [cm(1,:); [0.8 0.8 0.8]; cm(2,:); cm(4,:); cm(5,:)];

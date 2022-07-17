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
set(groot, 'defaultLineMarkerSize', 20);

odefile=@CO_oxidation;

parnames = {'k1', 'km1', 'k3', 'k2', 'km2', 'k4', 'lambda'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});
p(ind.k1) = 2.5;
p(ind.km1) = 1;
p(ind.k3) = 10;
p(ind.km2) = 0.1;
p(ind.k4) = 0.0675;
p(ind.lambda) = 1;
x  = [0.00295; 0.76211; 0.1678];

[x1_pred, v1_pred] = init_EP_EP(odefile, x, p, ind.k2);
opt = contset;
opt.MinStepsize   = 0.00001;
opt.InitStepsize  = 0.0001;
opt.MaxStepsize   = 0.1;
opt.Backward      = 1;
opt.MaxNumPoints  = 300;
opt.Singularities = 1;
[eqbr_x, ~, eqbr_bif_data] = cont(@equilibrium, x1_pred, v1_pred, opt);

%plot --width 1024 --height 800
hopfPointInfo   = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Hopf')==1);
limitPointsInfo = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Limit point')==1);
HopfPoint = eqbr_x(:,hopfPointInfo.index);
limitpoint1 = eqbr_x(:,limitPointsInfo(1).index);
limitpoint2 = eqbr_x(:,limitPointsInfo(2).index);
plot(eqbr_x(4,:), eqbr_x(3,:)); hold on
plot(HopfPoint(4), HopfPoint(3), '.r' ,'MarkerSize', 20)
plot(limitpoint1(4), limitpoint1(3),  '.k' ,'MarkerSize', 20)
plot(limitpoint2(4), limitpoint2(3),  '.k' ,'MarkerSize', 20)
xlabel('$k_2$')
ylabel('$s$')
legend({'Equilibrium curve', 'Hopf point', 'Limit points'}, 'Location', 'NorthWest')
title('Equilibrium curve in $(s,k_2)$-space')
axis([-16 2 0.05 0.35])

ap = [ind.k2 ind.lambda]; % continuation parameters
lp1.par = p;
lp1.par(ap) = eqbr_x(4, limitPointsInfo(1).index);
lp1.x = eqbr_x(1:3, limitPointsInfo(1).index);
[x1, v1] = init_LP_LP(odefile, lp1.x, lp1.par', ap);

opt.Backward = 0;
opt.TestTolerance = 1e-12;
opt.MaxTestIters = 10;
[lp_br1, ~, lp_br1_bif] = cont(@limitpoint, x1, v1, opt);

bt_points_info = lp_br1_bif(strcmp({lp_br1_bif.label}, 'BT')==1);

bt_index = bt_points_info(1).index;
bt1.x = lp_br1(1:3, bt_index);
bt1.par = p';
bt1.par(ap) = lp_br1(4:5, bt_index);
BToptions = BT_Hom_set_options();
BToptions.correct = false;
[x1_pred, v1_pred] = init_BT_Hom(odefile, bt1, ap, BToptions);

[hom1_x, hom1_v, ~] = newtcorr(x1_pred, v1_pred);
[x1_orbit, x1_saddle] = bt_rearr(hom1_x);

global homds
[homoclinic1_pred, saddle1_pred] = bt_rearr(x1_pred);
subplot(3,1,1); hold on;
title('Profiles of the predicted and correction homolinic orbits.')
plot(homds.finemsh, x1_orbit(1:3:end))
plot(homds.finemsh, homoclinic1_pred(1:3:end),'-.')
legend({'corrected', 'predicted'})
ylabel('$x$', 'fontsize', 12,'interpreter', 'latex')
subplot(3,1,2); hold on;
plot(homds.finemsh, x1_orbit(2:3:end))
plot(homds.finemsh, homoclinic1_pred(2:3:end),'-.')
legend({'corrected','predicted'})
ylabel('$y$', 'fontsize', 12,'interpreter', 'latex')
subplot(3,1,3); hold on;
plot(homds.finemsh, x1_orbit(3:3:end))
plot(homds.finemsh, homoclinic1_pred(3:3:end),'-.')
legend({'corrected',  'predicted'})
xlabel('$t$')
ylabel('$s$')

hold on
plot(x1_orbit(1:3:end),x1_orbit(2:3:end))
plot(homoclinic1_pred(1:3:end),homoclinic1_pred(2:3:end),'-.')
plot(x1_saddle(1), x1_saddle(2),'.', 'MarkerSize', 12, 'Color', [0 0.4470 0.7410])
plot(saddle1_pred(1), saddle1_pred(2),'.', 'MarkerSize', 12, 'Color', [0.8500, 0.3250, 0.0980])
xlabel('$x$')
ylabel('$y$')
title('Orbits and saddle points of predicted and corrected in phase-space')

opt.Singularities = 0;
opt.MinStepsize  = 1e-06; 
opt.InitStepsize = 1e-03; 
opt.MaxStepsize  = 1e-02; 
opt.MaxNumPoints = 300;
homoclinic_br1 = cont(@homoclinic, hom1_x, hom1_v, opt);

hold on
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
eps = linspace(0, 1.8);
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
plot(bt1.par(ind.k2), bt1.par(ind.lambda), '.k', 'MarkerSize', 20)
% set axis labels and legend
xlabel('$k_2$')
ylabel('$\lambda$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point'}, 'Location', 'NorthWest')
title('Comparision between computed and predicted parameter curve.')

hold on
%plot naive
plot(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br1(homds.coords(3:homds.nphase:end), 1:10:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$x$')
ylabel('$s$')
plot(bt1.x(1), bt1.x(3), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(x,s)$-phase space')

bt2_index = bt_points_info(2).index;
bt2.x = lp_br1(1:3, bt2_index);
bt2.par = p';
bt2.par(ap) = lp_br1(4:5, bt2_index);
options = BT_Hom_set_options();
[x2_pred, v2_pred] = init_BT_Hom(odefile, bt2,  ap, options);
homoclinic_br2 = cont(@homoclinic, x2_pred, v2_pred, opt);

[hopf_x, hopf_v] = init_BT_H(odefile, bt1.x, bt1.par, ap);

opt.MaxNumPoints = 2000;
opt.Singularities = 1;
opt.Eigenvalues = 1;
[hopf_br, ~, hopf_br_bif, ~, eigenvalues] = cont(@hopf, hopf_x, hopf_v, opt);
neutral_saddle_br = hopf_br(:,abs(imag(eigenvalues(2,:))) < 0.00001);
hopf_br_corrected = hopf_br(:,abs(imag(eigenvalues(2,:))) >= 0.00001);

hold on
colormap = lines();
homColor  = colormap(1,:);
hopfColor = colormap(2,:);
foldColor = colormap(5,:);
plot(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br1(homds.coords(3:homds.nphase:end), 1:10:end), ...
     'Color', homColor, 'HandleVisibility', 'Off')
plot(homoclinic_br2(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br2(homds.coords(3:homds.nphase:end), 1:10:end), ...
     'Color', homColor, 'HandleVisibility', 'Off')
plot(lp_br1(1,:), lp_br1(3,:), 'Color', foldColor)
plot(hopf_br_corrected(1,:), hopf_br_corrected(3,:), 'Color', hopfColor)
plot(neutral_saddle_br(1,1:end-1), neutral_saddle_br(3,1:end-1), '--', ...
    'Color', hopfColor)
xlabel('$x$')
ylabel('$s$')
plot(bt1.x(1), bt1.x(3), '.k')
plot(bt2.x(1), bt2.x(3), '.k', 'HandleVisibility', 'Off')
legend('Limit point branch', 'Hopf branch', 'Neutral saddle branch', ...
    'Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(x,s)$-phase space')
axis([0 0.17 0 0.7])

cusp_point_info = lp_br1_bif(strcmp({lp_br1_bif.label}, 'CP')==1);
genh_info = hopf_br_bif(strcmp({hopf_br_bif.label}, 'GH')==1);
cusp = lp_br1(:,cusp_point_info.index);
genh1 = hopf_br(:,genh_info(1).index);
genh2 = hopf_br(:,genh_info(2).index);

hold on
% plot computed parameter curve
plot(hopf_br_corrected(4,:), hopf_br_corrected(5,:), 'Color', hopfColor, 'linewidth', 1)
plot(lp_br1(4,:), lp_br1(5,:), 'Color', foldColor, 'linewidth', 1)
% plot homoclinic curves
plot(homoclinic_br1(homds.PeriodIdx+1,:), homoclinic_br1(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2)
plot(homoclinic_br2(homds.PeriodIdx+1,:), homoclinic_br2(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'HandleVisibility', 'Off', 'linewidth', 2)
% plot Bogdanov-Takens point
plot(bt1.par(ind.k2), bt1.par(ind.lambda), '.k')
plot(bt2.par(ind.k2), bt2.par(ind.lambda), '.k', 'HandleVisibility', 'Off')
plot(cusp(4), cusp(5), '.r')
plot(genh1(4), genh1(5), '.b')
plot(genh2(4), genh2(5), '.b', 'HandleVisibility', 'Off')
% set axis, labels and legend
xlabel('$k_2$')
ylabel('$\lambda$')
legend({'Limit point curve', 'Hopf curve', 'Homoclinic curve', ... 
    'Bogdanov-Takens points', 'Cusp point', 'Generalized Hopf point'}, ...
    'Location', 'NorthWest')
title('Comparision between computed and predicted parameter curve.')
axis([0.7130    1.4320    0.1413    1.0081])

set(groot, 'defaultLineMarkerSize', 7);
options = BT_Hom_set_options();
options.TTolerance = 1e-05;
options.messages = false;
options.correct = false;

amplitudes = logspace(-4, 0, 20);
methodList = {'orbital', 'LP', 'RegularPerturbation', ...
    'RegularPerturbationL2', 'LPHypernormalForm'};
relativeErrors = {};
for i=1:length(methodList)
    options.method = methodList{i};
    relativeErrors{i} = zeros(size(amplitudes));
    for j=1:length(amplitudes)
    options.amplitude = amplitudes(j);
    [x_pred, v0] = init_BT_Hom(odefile, bt1, ap, options);
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
title('CO-oxidation model')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
ax.ColorOrder = [cm(1,:); [0.8 0.8 0.8]; cm(2,:); cm(4,:); cm(5,:)];



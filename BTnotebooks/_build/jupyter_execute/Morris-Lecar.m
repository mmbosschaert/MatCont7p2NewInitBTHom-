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

odefile=@Morris_Lecar;

parnames = {'Iapp', 'v3'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});
p(ind.Iapp) = 0;
p(ind.v3) = 2;
x  = [-60.85456779; 0.0149139691];

[x1_pred, v1_pred] = init_EP_EP(odefile, x, p, ind.Iapp);
opt = contset;
opt.MaxStepsize   = 1;
opt.MaxNumPoints  = 300;
opt.Singularities = 1;
[eqbr_x, ~, eqbr_bif_data] = cont(@equilibrium, x1_pred, v1_pred, opt);

%plot --width 1024 --height 800
hopfPointsInfo   = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Hopf')==1);
HopfPoint1 = eqbr_x(:,hopfPointsInfo(1).index);
HopfPoint2 = eqbr_x(:,hopfPointsInfo(2).index);
plot(eqbr_x(3,:), eqbr_x(2,:)); hold on
plot(HopfPoint1(3), HopfPoint1(2), '.r' ,'MarkerSize', 20)
plot(HopfPoint2(3), HopfPoint2(2), '.r' ,'MarkerSize', 20)
xlabel('$I_{app}$')
ylabel('$w$')
legend({'Equilibrium curve', 'Hopf point'}, 'Location', 'NorthWest')
title('Equilibrium curve in $(I_{app},w)$-space')

ap = [ind.Iapp ind.v3]; % continuation parameters
hopf1.par = p;
hopf1.par(ap) = eqbr_x(3, hopfPointsInfo(1).index);
hopf1.x = eqbr_x(1:2, hopfPointsInfo(1).index);
[x1, v1] = init_H_H(odefile, hopf1.x, hopf1.par', ap);

opt.TestTolerance = 1e-15;
opt.MaxTestIters = 10;
opt.MaxNumPoints = 2000;
[hopf_br, ~, hopf_br_bif] = cont(@hopf, x1, v1, opt);

bt_points_info = hopf_br_bif(strcmp({hopf_br_bif.label}, 'BT')==1);
BTPoint1 = hopf_br(:,bt_points_info(1).index);
BTPoint2 = hopf_br(:,bt_points_info(2).index);
BTPoint3 = hopf_br(:,bt_points_info(3).index);
BTPoint4 = hopf_br(:,bt_points_info(4).index);
plot(hopf_br(3,:), hopf_br(4,:)); hold on
plot(BTPoint1(3), BTPoint1(4), '.b' ,'MarkerSize', 20)
plot(BTPoint2(3), BTPoint2(4), '.b' ,'MarkerSize', 20)
plot(BTPoint3(3), BTPoint3(4), '.b' ,'MarkerSize', 20)
plot(BTPoint4(3), BTPoint4(4), '.b' ,'MarkerSize', 20)
xlabel('$I_{app}$')
ylabel('$v_3$')
legend({'Hopf/Neutral saddle curve', 'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Hopf curve in $(I_{app},v_3)$-space')

bt_index = bt_points_info(1).index;
bt1.x = hopf_br(1:2, bt_index);
bt1.par = p';
bt1.par(ap) = hopf_br(3:4, bt_index);
BToptions = BT_Hom_set_options();
BToptions.correct = false;
[x1_pred, v1_pred] = init_BT_Hom(odefile, bt1, ap, BToptions);

[hom1_x, hom1_v, ~] = newtcorr(x1_pred, v1_pred);
[x1_orbit, x1_saddle] = bt_rearr(hom1_x);

[homoclinic1_pred, saddle1_pred] = bt_rearr(x1_pred);
subplot(2,1,1); hold on;
global homds
title('Profiles of the predicted and correction homolinic orbits.')
plot(homds.finemsh, x1_orbit(1:2:end))
plot(homds.finemsh, homoclinic1_pred(1:2:end),'-.')
legend({'corrected', 'predicted'})
ylabel('$V$')
subplot(2,1,2); hold on;
plot(homds.finemsh, x1_orbit(2:2:end))
plot(homds.finemsh, homoclinic1_pred(2:2:end),'-.')
legend({'corrected','predicted'})
ylabel('$w$')
legend({'corrected',  'predicted'})
xlabel('$t$')

hold on
plot(x1_orbit(1:2:end),x1_orbit(2:2:end))
plot(homoclinic1_pred(1:2:end),homoclinic1_pred(2:2:end),'-.')
plot(x1_saddle(1), x1_saddle(2),'.', 'MarkerSize', 12, 'Color', [0 0.4470 0.7410])
plot(saddle1_pred(1), saddle1_pred(2),'.', 'MarkerSize', 12, 'Color', [0.8500, 0.3250, 0.0980])
xlabel('$V$')
ylabel('$w$')
title('Orbits and saddle points of predicted and corrected in phase-space')

opt.Singularities = 0;
opt.MaxNumPoints = 1000;
homoclinic_br1 = cont(@homoclinic, hom1_x, hom1_v, opt);

hold on
% plot computed homoclinic parameter curve
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:));
% Bogdanov-Takens parameter-dependent smooth orbital normal form coefficients
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
plot(bt1.par(ind.Iapp), bt1.par(ind.v3), '.k', 'MarkerSize', 20)
% set axis labels and legend
xlabel('$I_{app}$')
ylabel('$v_3$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point'}, 'Location', 'SouthWest')
title('Comparision between computed and predicted parameter curve.')

hold on
plot(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br1(homds.coords(2:homds.nphase:end), 1:10:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
plot(bt1.x(1), bt1.x(2), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthWest')
title('Homoclic orbits in $(V,w)$-phase space')

options = BT_Hom_set_options();
options.messages = false;
options.correct = false;
options.TTolerance = 1.0e-05;

amplitudes = linspace(1.0e-03, 2,10);
XPredicted = zeros(330,length(amplitudes));
XCorrected = zeros(330,length(amplitudes));
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
plot(XPredicted(homds.coords(1:homds.nphase:end),1:10), ...
     XPredicted(homds.coords(2:homds.nphase:end),1:10), ...
      'color', cm(1,:))
plot(XCorrected(homds.coords(1:homds.nphase:end),1:10), ...
     XCorrected(homds.coords(2:homds.nphase:end),1:10), ...
      '--', 'color', cm(2,:))
plot(bt1.x(1), bt1.x(2), '.', 'MarkerSize', 16)
xlabel('$V$')
ylabel('$w$')
grid on

bt_index = bt_points_info(2).index;
bt2.x = hopf_br(1:2, bt_index);
bt2.par = p';
bt2.par(ap) = hopf_br(3:4, bt_index);
opt.MaxNumPoints = 100;
BToptions.correct = true;
[hom2_pred, hom2_v_pred] = init_BT_Hom(odefile, bt2, ap, BToptions);
homoclinic_br2 = cont(@homoclinic, hom2_pred, hom2_v_pred, opt);

hold on
plot(homoclinic_br2(homds.coords(1:homds.nphase:end), 1:5:end), ...
    homoclinic_br2(homds.coords(2:homds.nphase:end), 1:5:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
plot(bt2.x(1), bt2.x(2), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthWest')
title('Homoclic orbits in $(V,w)$-phase space')

bt_index = bt_points_info(3).index;
bt3.x = hopf_br(1:2, bt_index);
bt3.par = p';
bt3.par(ap) = hopf_br(3:4, bt_index);
[hom3_pred, hom3_v_pred] = init_BT_Hom(odefile, bt3, ap, BToptions);
opt.MaxNumPoints  = 200;
homoclinic_br3 = cont(@homoclinic, hom3_pred, hom3_v_pred, opt);

hold on
plot(homoclinic_br3(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br3(homds.coords(2:homds.nphase:end), 1:10:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
plot(bt3.x(1), bt3.x(2), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthWest')
title('Homoclic orbits in $(V,w)$-phase space')

bt_index = bt_points_info(4).index;
bt4.x = hopf_br(1:2, bt_index);
bt4.par = p';
bt4.par(ap) = hopf_br(3:4, bt_index);
[hom4_pred, hom4_v_pred] = init_BT_Hom(odefile, bt4, ap, BToptions);
homoclinic_br4 = cont(@homoclinic, hom4_pred, hom4_v_pred, opt);

hold on
plot(homoclinic_br4(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br4(homds.coords(2:homds.nphase:end), 1:10:end), ...
     'Color', [0 0.4470 0.7410], 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
plot(bt4.x(1), bt4.x(2), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthWest')
title('Homoclic orbits in $(V,w)$-phase space')

[lp1_x, lp1_v] = init_BT_LP(odefile, bt1.x, bt1.par, ap);
opt.MaxNumPoints = 2000;
opt.Singularities = 1;
opt.Backward = 1;
[lp1_br, ~, lp1_br_bif] = cont(@limitpoint, lp1_x, lp1_v, opt);
opt.Backward = 0;
opt.MaxNumPoints = 200;
lp1_br_rev = cont(@limitpoint, lp1_x, lp1_v, opt);

[lp2_x, lp2_v] = init_BT_LP(odefile, bt3.x, bt3.par, ap);
opt.MaxNumPoints = 2000;
opt.Singularities = 1;
opt.Backward = 0;
[lp2_br, ~, lp2_br_bif] = cont(@limitpoint, lp2_x, lp2_v, opt);
opt.Backward = 1;
opt.MaxNumPoints = 200;
lp2_br_rev = cont(@limitpoint, lp2_x, lp2_v, opt);

cusp1_info = lp1_br_bif(strcmp({lp1_br_bif.label}, 'CP')==1);
cusp2_info = lp2_br_bif(strcmp({lp2_br_bif.label}, 'CP')==1);
Cusp1 = lp1_br(:,cusp1_info.index);
Cusp2 = lp2_br(:,cusp2_info.index);

%plot --width 1024 --height 800
figure; hold on 
homColor  = cm(1,:);
hopfColor = cm(2,:);
foldColor = cm(5,:);
plot(hopf_br(1,:), hopf_br(2,:), 'Color', hopfColor); hold on
plot(lp2_br(1,:), lp2_br(2,:), 'Color', foldColor);
plot(lp2_br_rev(1,:), lp2_br_rev(2,:), 'Color', foldColor, ...
    'HandleVisibility', 'Off');
plot(Cusp2(1), Cusp2(2), '.r' ,'MarkerSize', 20)
plot(BTPoint2(1), BTPoint2(2), '.k' ,'MarkerSize', 20)
plot(BTPoint3(1), BTPoint3(2), '.k' ,'MarkerSize', 20)
plot(homoclinic_br2(homds.coords(1:homds.nphase:end), 1:5:end), ...
     homoclinic_br2(homds.coords(2:homds.nphase:end), 1:5:end), ...
     'Color', homColor, 'DisplayName', 'Off')
plot(homoclinic_br3(homds.coords(1:homds.nphase:end), 1:5:end), ...
     homoclinic_br3(homds.coords(2:homds.nphase:end), 1:5:end), ...
     'Color', homColor, 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
legend({'Hopf/Neutral saddle curve', 'Limit point curve', 'Cusp point', ...
    'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Hopf curve in $(V,w)$-space')
axis([-21.7891   12.0474    0.8746    1.0000])

plot(hopf_br(1,:), hopf_br(2,:), 'Color', hopfColor); hold on
plot(lp1_br(1,:), lp1_br(2,:), 'Color', foldColor);
plot(lp1_br_rev(1,:), lp1_br_rev(2,:), 'Color', foldColor, ...
    'HandleVisibility', 'Off');
plot(Cusp1(1), Cusp1(2), '.r' ,'MarkerSize', 20)
plot(BTPoint1(1), BTPoint1(2), '.k' ,'MarkerSize', 20)
plot(BTPoint4(1), BTPoint4(2), '.k' ,'MarkerSize', 20)
plot(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br1(homds.coords(2:homds.nphase:end), 1:10:end), ...
     'Color', homColor, 'HandleVisibility', 'Off')
plot(homoclinic_br4(homds.coords(1:homds.nphase:end), 1:10:end), ...
     homoclinic_br4(homds.coords(2:homds.nphase:end), 1:10:end), ...
     'Color', homColor, 'HandleVisibility', 'Off')
xlabel('$V$')
ylabel('$w$')
legend({'Hopf/Neutral saddle curve', 'Limit point curve', 'Cusp point', ...
    'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Hopf curve in $(V,w)$-space')
axis([-50 50 -0.05 0.45])

hold on
plot(hopf_br(3,:), hopf_br(4,:), 'Color', hopfColor, 'linewidth', 2)
plot(lp1_br(3,:), lp1_br(4,:), 'Color', foldColor, 'linewidth', 2)
plot(lp1_br_rev(3,:), lp1_br_rev(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(lp2_br(3,:), lp2_br(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(lp2_br_rev(3,:), lp2_br_rev(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br2(homds.PeriodIdx+1,:), ...
     homoclinic_br2(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br3(homds.PeriodIdx+1,:), ...
     homoclinic_br3(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br4(homds.PeriodIdx+1,:), ...
     homoclinic_br4(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(Cusp2(3), Cusp2(4), '.r' ,'MarkerSize', 20)
plot(BTPoint1(3), BTPoint1(4), '.b' ,'MarkerSize', 20)
plot(BTPoint2(3), BTPoint2(4), '.b' ,'MarkerSize', 20)
plot(BTPoint3(3), BTPoint3(4), '.b' ,'MarkerSize', 20)
plot(BTPoint4(3), BTPoint4(4), '.b' ,'MarkerSize', 20)
xlabel('$I_{app}$')
ylabel('$v_3$')
legend({'Hopf/Neutral Saddle curve', 'Fold curve', 'Cusp point',...
    'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Bifurcation daigram in $(I_{app},v_3)$-space')
axis([456.5097  519.9544  -63.9705  -34.7564])

hold on
plot(hopf_br(3,:), hopf_br(4,:), 'Color', hopfColor, 'linewidth', 2)
plot(lp1_br(3,:), lp1_br(4,:), 'Color', foldColor, 'linewidth', 2)
plot(lp1_br_rev(3,:), lp1_br_rev(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(lp2_br(3,:), lp2_br(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(lp2_br_rev(3,:), lp2_br_rev(4,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off')
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br2(homds.PeriodIdx+1,:), ...
     homoclinic_br2(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br3(homds.PeriodIdx+1,:), ...
     homoclinic_br3(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(homoclinic_br4(homds.PeriodIdx+1,:), ...
     homoclinic_br4(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
plot(Cusp1(3), Cusp1(4), '.r' ,'MarkerSize', 20)
plot(BTPoint1(3), BTPoint1(4), '.b' ,'MarkerSize', 20)
plot(BTPoint2(3), BTPoint2(4), '.b' ,'MarkerSize', 20)
plot(BTPoint3(3), BTPoint3(4), '.b' ,'MarkerSize', 20)
plot(BTPoint4(3), BTPoint4(4), '.b' ,'MarkerSize', 20)
xlabel('$I_{app}$')
ylabel('$v_3$')
legend({'Hopf/Neutral Saddle curve', 'Fold curve', 'Cusp point',...
    'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Bifurcation daigram in $(I_{app},v_3)$-space')
axis([-231.3697  104.5793    4.0404   73.0235])

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

loglog(amplitudes, relativeErrors{1}(:), 'd', ...
       amplitudes, relativeErrors{2}(:), '--', ...
       amplitudes, relativeErrors{3}(:), '*', ...
       amplitudes, relativeErrors{4}(:), 's', ...
       amplitudes, relativeErrors{5}(:), '+')
legend(methodList, 'Location', 'NorthWest')
title('Morris-Lecar model')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
ax.ColorOrder = [cm(1,:); [0.8 0.8 0.8]; cm(2,:); cm(4,:); cm(5,:)];

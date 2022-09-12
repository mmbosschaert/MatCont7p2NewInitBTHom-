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

odefile=@HomoclinicRGflows;

parnames = {'epsilon', 'M', 'N'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});
p(ind.epsilon) = 1;
p(ind.M) = 0.2945;
p(ind.N) = 4.036;
x  = [0.0701457361241472, -0.06520883770451065, 0.001823543197553845, 0.22874527306411319]';

[x1_pred, v1_pred] = init_EP_EP(odefile, x, p, ind.M);
opt = contset;
opt.MaxNumPoints  = 300;
opt.Singularities = 1;
opt.Backward = 1;
[eqbr_x, ~, eqbr_bif_data] = cont(@equilibrium, x1_pred, v1_pred, opt);

%plot --width 1024 --height 800
plot(eqbr_x(5,:), eqbr_x(4,:)); hold on
foldInfo   = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Limit point')==1);
foldInfocell = struct2cell(foldInfo);
foldInd = cell2mat(foldInfocell(1,:));
hopfInfo   = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Hopf')==1);
hopfInfocell = struct2cell(hopfInfo);
hopfInd = cell2mat(hopfInfocell(1,:));
plot(eqbr_x(5,foldInd), eqbr_x(4,foldInd), '.r', 'MarkerSize', 20); hold on
plot(eqbr_x(5,hopfInd), eqbr_x(4,hopfInd), '.b', 'MarkerSize', 20); hold on
xlabel('$M$')
ylabel('$g_4$')
legend({'Equilibrium curve'}, 'Location', 'NorthEast')
title('Equilibrium curve in $(M,g_4)$-space')

ap = [ind.M ind.N];
hopfInfo = eqbr_bif_data(strcmp({eqbr_bif_data.msg}, 'Hopf')==1);
hopf.x = eqbr_x(1:4,hopfInfo(2).index);
hopf.par = p';
hopf.par(ind.M) = eqbr_x(5,hopfInfo(2).index);
[hopf1_x, hopf1_v] = init_H_H(odefile, hopf.x, hopf.par, ap);

opt.TestTolerance = 1e-12;
opt.MaxTestIters = 10;
opt.Backward = 0;
opt.MaxNumPoints = 50;
[hopf_br, ~, hopf_br_bif] = cont(@hopf, hopf1_x, hopf1_v, opt);

bt_points_info = hopf_br_bif(strcmp({hopf_br_bif.label}, 'BT')==1);
BTPoint1 = hopf_br(:,bt_points_info(1).index);
BTPoint2 = hopf_br(:,bt_points_info(2).index);
plot(hopf_br(5,:), hopf_br(6,:)); hold on
plot(BTPoint1(5), BTPoint1(6), '.b' ,'MarkerSize', 20)
plot(BTPoint2(5), BTPoint2(6), '.b' ,'MarkerSize', 20)
xlabel('$M$')
ylabel('$N$')
legend({'Hopf branch', 'Bogadanov-Takens point'}, 'Location', 'NorthWest')
title('Hopf curve in $(M,N)$-space')

bt_index = bt_points_info(1).index;
bt1.x = hopf_br(1:4, bt_index);
bt1.par = p';
bt1.par(ap) = hopf_br(5:6, bt_index);
BToptions = BT_Hom_set_options();
BToptions.correct = false;
BToptions.amplitude = 0.2;
[x1_pred, v1_pred] = init_BT_Hom(odefile, bt1, ap, BToptions);

[hom1_x, hom1_v, ~] = newtcorr(x1_pred, v1_pred);
[x1_orbit, x1_saddle] = bt_rearr(hom1_x);

[homoclinic1_pred, saddle1_pred] = bt_rearr(x1_pred);
subplot(4,1,1); hold on;
global homds
title('Profiles of the predicted and correction homolinic orbits.')
plot(homds.finemsh, x1_orbit(1:4:end))
plot(homds.finemsh, homoclinic1_pred(1:4:end),'.')
legend({'corrected', 'predicted'})
ylabel('$g_1$')
subplot(4,1,2); hold on;
plot(homds.finemsh, x1_orbit(2:4:end))
plot(homds.finemsh, homoclinic1_pred(2:4:end),'.')
legend({'corrected','predicted'})
ylabel('$g_2$')
subplot(4,1,3); hold on;
plot(homds.finemsh, x1_orbit(3:4:end))
plot(homds.finemsh, homoclinic1_pred(3:4:end),'.')
legend({'corrected','predicted'})
ylabel('$g_3$')
subplot(4,1,4); hold on;
plot(homds.finemsh, x1_orbit(4:4:end))
plot(homds.finemsh, homoclinic1_pred(4:4:end),'.')
legend({'corrected','predicted'})
ylabel('$g_4$')
legend({'corrected',  'predicted'})
xlabel('$t$')

hold on
plot(x1_orbit(1:4:end),x1_orbit(2:4:end))
plot(homoclinic1_pred(1:4:end),homoclinic1_pred(2:4:end),'.')
plot(x1_saddle(1), x1_saddle(2),'.', 'MarkerSize', 12, 'Color', [0 0.4470 0.7410])
plot(saddle1_pred(1), saddle1_pred(2),'.', 'MarkerSize', 12, 'Color', [0.8500, 0.3250, 0.0980])
xlabel('$g_1$')
ylabel('$g_2$')
title('Orbits and saddle points of predicted and corrected in phase-space')

[homoclinic_br1, homoclinic_br1_v, homoclinic_singularities] = cont(@homoclinic, hom1_x, hom1_v, opt);

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
eps = linspace(0, 0.05);
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
plot(bt1.par(ind.M), bt1.par(ind.N), '.k', 'MarkerSize', 20)
% set axis labels and legend
xlabel('$M$')
ylabel('$N$')
legend({'Homoclinic curve', 'Current homoclinic predictor', ...
    'Bogdanov-Takens point'}, 'Location', 'SouthEast')
title('Comparision between computed and predicted parameter curve.')

global homds
cm = lines;
hold on
plot3(homoclinic_br1(homds.coords(1:homds.nphase:end), 1:4:end), ...
      homoclinic_br1(homds.coords(2:homds.nphase:end), 1:4:end), ...
      homoclinic_br1(homds.coords(3:homds.nphase:end), 1:4:end), ...
      'Color', cm(1,:), 'HandleVisibility', 'Off')
bif_points = struct2cell(homoclinic_singularities);
plot3(homoclinic_br1(homds.coords(1:homds.nphase:end), cell2mat(bif_points(1,2:end-1))), ...
      homoclinic_br1(homds.coords(2:homds.nphase:end), cell2mat(bif_points(1,2:end-1))), ...
      homoclinic_br1(homds.coords(3:homds.nphase:end), cell2mat(bif_points(1,2:end-1))), ...
      'Color', cm(2,:), 'HandleVisibility', 'Off', 'LineWidth', 2)
xlabel('$g_1$')
ylabel('$g_2$')
zlabel('$g_3$')
plot3(bt1.x(1), bt1.x(2), bt1.x(3), '.k' ,'MarkerSize', 20)
legend('Bogdanov-Takens point', 'Location', 'SouthEast')
title('Homoclic orbits in $(g_1,g_2,g_3)$-phase space')
grid on
view(66, 50) 

options = BT_Hom_set_options();
options.messages = false;
options.correct = false;
options.TTolerance = 1.0e-05;

amplitudes = linspace(1.0e-03, 1.0e-01, 10);
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

clf
subplot(2,2,1); hold on
R = @(alpha) [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
S = @(s) [s 0; 0 1];
parsCorrected = XCorrected(homds.PeriodIdx+1,1:end).*ones(homds.tps,10);
for i=1:length(amplitudes)
    [profile, saddle] = bt_rearr(XCorrected(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(100)*R(1.2181)*(profile(1:2,:) ...
                        - saddle(1:2)) + saddle(1:2));
    plot3(parsCorrected(:,i),profileRotated(1,:)', profileRotated(2,:)', ...
          '.','color', cm(1,:))

    [profile, saddle] = bt_rearr(XPredicted(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(100)*R(1.2181)*(profile(1:2,:) ...
                        - saddle(1:2)) + saddle(1:2));
    plot3(parsCorrected(:,i),profileRotated(1,:)', profileRotated(2,:)', ...
          'color', cm(2,:))
end
xlabel('$M$')
ylabel('$\tilde g_1$')
zlabel('$\tilde g_2$')
grid on
view(32,15)

subplot(2,2,2); hold on
for i=1:length(amplitudes)
    alpha0 = -1.139;
    [profile, saddle] = bt_rearr(XCorrected(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-alpha0)*(S(100)*R(alpha0)*(profile(3:4,:) ...
                        - saddle(3:4)) + saddle(3:4));
    plot3(parsCorrected(:,i),profileRotated(1,:)', profileRotated(2,:)', ...
          '.','color', cm(1,:))

    [profile, saddle] = bt_rearr(XPredicted(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-alpha0)*(S(100)*R(alpha0)*(profile(3:4,:) ...
                        - saddle(3:4)) + saddle(3:4));
    plot3(parsCorrected(:,i),profileRotated(1,:)', profileRotated(2,:)', ...
          'color', cm(2,:))
end
xlabel('$M$')
ylabel('$\tilde g_3$')
zlabel('$\tilde g_4$')
grid on
view(28,15)

subplot(2,2,3); hold on
R = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
S = @(s) [s 0 0; 0 1 0; 0 0 1];
for i=1:length(amplitudes)
    [profile, saddle] = bt_rearr(XCorrected(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(200)*R(1.2181)*(profile(1:3,:) ...
                        - saddle(1:3)) + saddle(1:3));
    plot3(profileRotated(1,:)', profileRotated(2,:)', profileRotated(3,:)', ...
          '.','color', cm(1,:))

    [profile, saddle] = bt_rearr(XPredicted(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(200)*R(1.2181)*(profile(1:3,:) ...
                        - saddle(1:3)) + saddle(1:3));
    plot3(profileRotated(1,:)', profileRotated(2,:)', profileRotated(3,:)', ...
          'color', cm(2,:))
end
xlabel('$\tilde g_1$')
ylabel('$\tilde g_2$')
zlabel('$g_3$')
grid on
view(332,11)

subplot(2,2,4); hold on
for i=1:length(amplitudes)
    [profile, saddle] = bt_rearr(XCorrected(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(200)*R(1.2181)*(profile([1,2,4],:) ...
                        - saddle([1,2,4])) + saddle([1,2,4]));
    plot3(profileRotated(1,:)', profileRotated(2,:)', profileRotated(3,:)', ...
          '.','color', cm(1,:))

    [profile, saddle] = bt_rearr(XPredicted(:,i));
    profile = reshape(profile,4,[]);
    profileRotated = R(-1.2181)*(S(200)*R(1.2181)*(profile([1,2,4],:) ...
                        - saddle([1,2,4])) + saddle([1,2,4]));
    plot3(profileRotated(1,:)', profileRotated(2,:)', profileRotated(3,:)', ...
          'color', cm(2,:))
end
xlabel('$\tilde g_1$')
ylabel('$\tilde g_2$')
zlabel('$g_4$')
grid on
view(332,11)

[lp1_x, lp1_v] = init_BT_LP(odefile, bt1.x, bt1.par, ap);
[lp_br, ~, lp_br1_bif] = cont(@limitpoint, lp1_x, lp1_v, opt);
opt.Backward = 1;
lp_br_rev = cont(@limitpoint, lp1_x, lp1_v, opt);

%plot inline 
hold on
homColor  = cm(1,:);
hopfColor = cm(2,:);
foldColor = cm(5,:);
plot(hopf_br(5,:), hopf_br(6,:), 'Color', hopfColor, 'linewidth', 2)
plot(lp_br(5,:), lp_br(6,:), 'Color', foldColor, 'linewidth', 2)
plot(lp_br_rev(5,:), lp_br_rev(6,:), 'Color', foldColor, 'linewidth', 2, ...
    'HandleVisibility', 'Off', 'linewidth', 2);
plot(BTPoint1(5), BTPoint1(6), '.b' ,'MarkerSize', 20)
plot(homoclinic_br1(homds.PeriodIdx+1,:), ...
     homoclinic_br1(homds.PeriodIdx+2,:), ...
     '--', 'Color', homColor, 'linewidth', 2, 'HandleVisibility', 'Off')
xlabel('$M$')
ylabel('$N$')
legend({'Hopf/Neutral Saddle curve', 'Fold curve',...
    'Bogadanov-Takens point'}, 'Location', 'NorthEast')
title('Bifurcation daigram in $(M,N)$-space')
axis([0.1956    0.4852    3.9896    4.1020])

BToptions = BT_Hom_set_options();
BToptions.TTolerance = 1e-05;
BToptions.messages = false;
BToptions.correct = false;

amplitudes = logspace(-4, 0, 20);
methodList = {'orbital', 'LP', 'RegularPerturbation', ...
    'RegularPerturbationL2', 'LPHypernormalForm'};
relativeErrors = {};
for i=1:length(methodList)
    for o=1:3
        BToptions.method = methodList{i};
        BToptions.order = o;
        relativeErrors{o,i} = zeros(size(amplitudes));
        for j=1:length(amplitudes)
            BToptions.amplitude = amplitudes(j);
            [x_pred, v0] = init_BT_Hom(odefile, bt1, ap, BToptions);
            try
                x_corrected = newtcorr(x_pred, v0);
                relativeErrors{o,i}(j) = norm(x_corrected-x_pred)/norm(x_corrected);
            catch
                warning('Did not converge.')
                continue
            end
        end
    end
end

cm = lines();
loglog(amplitudes, relativeErrors{3,1}(:), 'd', ...
       amplitudes, relativeErrors{3,2}(:), '--', ...
       amplitudes, relativeErrors{3,3}(:), '*', ...
       amplitudes, relativeErrors{3,4}(:), 's', ...
       amplitudes, relativeErrors{3,5}(:), '+')
legend(methodList, 'Location', 'NorthWest')
title('Hodgkin-Huxley equations')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
ax.ColorOrder = [cm(1,:); [0.8 0.8 0.8]; cm(2,:); cm(4,:); cm(5,:)];

writematrix([amplitudes', relativeErrors{1,3}(:)], '../../data/HomRGflowsRPorder1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{2,3}(:)], '../../data/HomRGflowsRPorder2.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{3,3}(:)], '../../data/HomRGflowsRPorder3.csv', 'Delimiter', ' ')
                                                   
writematrix([amplitudes', relativeErrors{1,2}(:)], '../../data/HomRGflowsLPorder1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{2,2}(:)], '../../data/HomRGflowsLPorder2.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{3,2}(:)], '../../data/HomRGflowsLPorder3.csv', 'Delimiter', ' ')

writematrix([amplitudes', relativeErrors{3,1}(:)], '../../data/HomRGflowsLPorder3orbital.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{3,4}(:)], '../../data/HomRGflowsRegularPerturbationL2.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors{3,5}(:)], '../../data/HomRGflowsLPHypernormalForm.csv', 'Delimiter', ' ')

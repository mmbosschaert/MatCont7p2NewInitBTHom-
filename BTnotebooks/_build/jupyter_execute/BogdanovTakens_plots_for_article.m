clear all
matcontpath = '../';
addpath(matcontpath)
addpath([matcontpath, 'Systems'])
addpath([matcontpath, 'Equilibrium'])
addpath([matcontpath, 'LimitPoint'])
addpath([matcontpath, 'LimitPointCycle'])
addpath([matcontpath, 'Hopf'])
addpath([matcontpath, 'Homoclinic'])
addpath([matcontpath, 'LimitCycle'])
addpath([matcontpath, 'Continuer'])
addpath([matcontpath, 'MultilinearForms'])
addpath([matcontpath, 'Utilities'])
set(groot, 'defaultTextInterpreter', 'LaTeX');
set(0,'defaultAxesFontSize',15)

odefile=@BogdanovTakensNormalForm;

w0 = 0;
w1 = 0;
beta1 = 0;
beta2 = 0;
bt.x = [w0; w1];
bt.par = [beta1; beta2];

parnames = {'beta1', 'beta2'};
cind = [parnames;num2cell(1:length(parnames))];
ind  = struct(cind{:});

%plot --width 1024 --height 800
ap = [ind.beta1, ind.beta2];
BToptions = BT_Hom_set_options();
BToptions.method = 'LP';
BToptions.HigherOrderTimeReparametrization = 0;
BToptions.correct = 0;
BToptions.amplitude = 6;
BToptions.TTolerance = 1.0e-03;
hom2016pred = init_BT_Hom(odefile, bt,  ap, BToptions);
homcorrected = newtcorr(hom2016pred, hom_v_pred);
global homds
[hom2016predOrbit, saddle2016pred] = bt_rearr(hom2016pred);
[homcorrectedOrbit, saddlecorrected] = bt_rearr(homcorrected);
subplot(2,1,1); hold on;
title('Profiles of the predicted and correction homolinic orbits.')
finemsh = 6*(2*homds.finemsh-1);
plot(finemsh, homcorrectedOrbit(1:2:end),'-')
plot(finemsh, hom2016predOrbit(1:2:end),'--')
legend({'corrected', 'predicted 2016'})
ylabel('$w_0$')
subplot(2,1,2); hold on;
plot(finemsh, homcorrectedOrbit(2:2:end),'-')
plot(finemsh, hom2016predOrbit(2:2:end),'--')
legend({'corrected', 'predicted 2016'})
ylabel('$w_1$')
xlabel('$t$')

% write data to file
writematrix([finemsh', homcorrectedOrbit(1:2:end)], '../../data/w0CorrectedLPWithoutTimeInversion3.csv', 'Delimiter', ' ')
writematrix([finemsh', hom2016predOrbit(1:2:end)], '../../data/w0PredictorLPWithoutTimeInversion3.csv', 'Delimiter', ' ')
writematrix([finemsh', homcorrectedOrbit(2:2:end)], '../../data/w1CorrectedLPWithoutTimeInversion3.csv', 'Delimiter', ' ')
writematrix([finemsh', hom2016predOrbit(2:2:end)], '../../data/w1PredictorLPWithoutTimeInversion3.csv', 'Delimiter', ' ')

%plot --width 1024 --height 800
ap = [ind.beta1, ind.beta2];
BToptions = BT_Hom_set_options();
BToptions.method = 'RegularPerturbation';
BToptions.correct = 0;
BToptions.amplitude = 6;
BToptions.TTolerance = 1.0e-03;
[homRPpred, hom_v_pred] = init_BT_Hom(odefile, bt,  ap, BToptions);
homRPcorrected = newtcorr(homRPpred, hom_v_pred);
global homds
[homRPpredOrbit, ~] = bt_rearr(homRPpred);
[homRPcorrectedOrbit, ~] = bt_rearr(homRPcorrected);
subplot(2,1,1); hold on;
title('Profiles of the predicted and correction homolinic orbits.')
finemsh = 6*(2*homds.finemsh-1);
plot(finemsh, homRPcorrectedOrbit(1:2:end),'-')
plot(finemsh, homRPpredOrbit(1:2:end),'--')
legend({'corrected', 'predicted RP'})
ylabel('$w_0$')
subplot(2,1,2); hold on;
plot(finemsh, homRPcorrectedOrbit(2:2:end),'-')
plot(finemsh, homRPpredOrbit(2:2:end),'--')
legend({'corrected', 'predicted RP'})
ylabel('$w_1$')
xlabel('$t$')

% write data to file
writematrix([finemsh', homRPcorrectedOrbit(1:2:end)], '../../data/w0CorrectedRP3.csv', 'Delimiter', ' ')
writematrix([finemsh', homRPpredOrbit(1:2:end)], '../../data/w0PredictorRP3.csv', 'Delimiter', ' ')
writematrix([finemsh', homRPcorrectedOrbit(2:2:end)], '../../data/w1CorrectedRP3.csv', 'Delimiter', ' ')
writematrix([finemsh', homRPpredOrbit(2:2:end)], '../../data/w1PredictorRP3.csv', 'Delimiter', ' ')

%plot --width 1024 --height 800
ap = [ind.beta1, ind.beta2];
BToptions = BT_Hom_set_options();
BToptions.method = 'LP';
BToptions.correct = 0;
BToptions.amplitude = 6;
BToptions.TTolerance = 1.0e-03;
[hom2021pred, hom_v_pred] = init_BT_Hom(odefile, bt,  ap, BToptions);
homcorrected = newtcorr(hom2021pred, hom_v_pred);
global homds
[hom2021predOrbit, saddle2021pred] = bt_rearr(hom2021pred);
[homcorrectedOrbit, saddlecorrected] = bt_rearr(homcorrected);
subplot(2,1,1); hold on;
title('Profiles of the predicted and correction homolinic orbits.')
finemsh = 6*(2*homds.finemsh-1);
plot(finemsh, homcorrectedOrbit(1:2:end),'-')
plot(finemsh, hom2021predOrbit(1:2:end),'--')
legend({'corrected', 'predicted 2021'})
ylabel('$w_0$')
subplot(2,1,2); hold on;
plot(finemsh, homcorrectedOrbit(2:2:end),'-')
plot(finemsh, hom2021predOrbit(2:2:end),'--')
legend({'corrected', 'predicted 2021'})
ylabel('$w_1$')
xlabel('$t$')

% write data to file
writematrix([finemsh', homcorrectedOrbit(1:2:end)], '../../data/w0CorrectedLP3.csv', 'Delimiter', ' ')
writematrix([finemsh', hom2021predOrbit(1:2:end)], '../../data/w0PredictorLP3.csv', 'Delimiter', ' ')
writematrix([finemsh', homcorrectedOrbit(2:2:end)], '../../data/w1CorrectedLP3.csv', 'Delimiter', ' ')
writematrix([finemsh', hom2021predOrbit(2:2:end)], '../../data/w1PredictorLP3.csv', 'Delimiter', ' ')

hold on
plot(homcorrectedOrbit(1:2:end),homcorrectedOrbit(2:2:end), '-')
plot(hom2016predOrbit(1:2:end),hom2016predOrbit(2:2:end),'--')
plot(hom2021predOrbit(1:2:end),hom2021predOrbit(2:2:end),'x')
plot(saddlecorrected(1), saddlecorrected(2), '.', 'MarkerSize', 12, 'Color', [0 0.4470 0.7410])
xlabel('$w_0$')
ylabel('$w_1$')
legend({'corrected', 'predicted 2016', 'predicted 2021', 'Corrected saddle point'})
title('Orbits and saddle points of predicted and corrected in phase-space')

BToptions = BT_Hom_set_options();
BToptions.TTolerance = 1e-05;
BToptions.messages = false;
BToptions.correct = false;

amplitudes = logspace(-3.5, 0, 40);
methodList = {'RegularPerturbation', 'RegularPerturbationL2', 'LP', 'orbitalv2'};

relativeErrors = {};
relativeErrors_w0 = {};
relativeErrors_w1 = {};
for i=1:4
    for k=1:length(methodList)
        relativeErrors{i,k} = zeros(size(amplitudes));
        relativeErrors_w0{i,k} = zeros(size(amplitudes));
        relativeErrors_w1{i,k} = zeros(size(amplitudes));
        BToptions.order = i-1;
        BToptions.method = methodList{k}; 
        for j=1:length(amplitudes)
            BToptions.amplitude = amplitudes(j);
            [x_pred, v0] = init_BT_Hom(odefile, bt, ap, BToptions);
            try
                x_corrected = newtcorr(x_pred, v0);
                relativeErrors{i,k}(j) = norm(x_corrected-x_pred)/norm(x_corrected);
                [hom_corrected_orbit, ~] = bt_rearr(x_corrected);
                [hompred_orbit, ~] = bt_rearr(x_pred);
                relativeErrors_w0{i,k}(j) = norm(hom_corrected_orbit(1:2:end) - hompred_orbit(1:2:end))/norm(hom_corrected_orbit(1:2:end));
                relativeErrors_w1{i,k}(j) = norm(hom_corrected_orbit(2:2:end) - hompred_orbit(2:2:end))/norm(hom_corrected_orbit(2:2:end));
            catch
                warning('Did not converge.')
                continue
            end
        end
    end
end

i=1;
k=5;
relativeErrors{i} = zeros(size(amplitudes));
BToptions.HigherOrderTimeReparametrization = 0;
BToptions.method = 'LP';
for j=1:length(amplitudes)
    BToptions.amplitude = amplitudes(j);
    [x_pred, v0] = init_BT_Hom(odefile, bt, ap, BToptions);
    try
        x_corrected = newtcorr(x_pred, v0);
        relativeErrors{i}(j) = norm(x_corrected-x_pred)/norm(x_corrected);
        [hom_corrected_orbit, ~] = bt_rearr(x_corrected);
        [hompred_orbit, ~] = bt_rearr(x_pred);
        relativeErrors_w0{i,k}(j) = norm(hom_corrected_orbit(1:2:end) - hompred_orbit(1:2:end))/norm(hom_corrected_orbit(1:2:end));
        relativeErrors_w1{i,k}(j) = norm(hom_corrected_orbit(2:2:end) - hompred_orbit(2:2:end))/norm(hom_corrected_orbit(2:2:end));
    catch
        warning('Did not converge.')
        continue
    end
end

cm = lines();
loglog(amplitudes, relativeErrors_w0{1,1}(:), 'o', ...
       amplitudes, relativeErrors_w0{2,1}(:), '+', ...
       amplitudes, relativeErrors_w0{3,1}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,1}(:), '*', ...
       amplitudes, relativeErrors_w0{2,2}(:), '+', ...
       amplitudes, relativeErrors_w0{3,2}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,2}(:), '*')
legend('Regular Perturbation order 0','Regular Perturbation order 1','Regular Perturbation order 2', ...
'Regular Perturbation order 3', ...
'Regular Perturbation L2 phase condition order 1', ...
'Regular Perturbation L2 phase condition order 2', ...
'Regular Perturbation L2 phase condition order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [cm(1,:); cm(1,:); cm(1,:); cm(1,:); ...
                 cm(2,:); cm(2,:); cm(2,:)];

loglog(amplitudes, relativeErrors_w1{1,1}(:), 'o', ...
       amplitudes, relativeErrors_w1{2,1}(:), '+', ...
       amplitudes, relativeErrors_w1{3,1}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,1}(:), '*', ...
       amplitudes, relativeErrors_w1{2,2}(:), '+', ...
       amplitudes, relativeErrors_w1{3,2}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,2}(:), '*')
legend('Regular Perturbation order 0','Regular Perturbation order 1','Regular Perturbation order 2', ...
'Regular Perturbation order 3', ...
'Regular Perturbation L2 phase condition order 1', ...
'Regular Perturbation L2 phase condition order 2', ...
'Regular Perturbation L2 phase condition order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [cm(1,:); cm(1,:); cm(1,:); cm(1,:); ...
                 cm(2,:); cm(2,:); cm(2,:)];

loglog(amplitudes, relativeErrors_w0{1,1}(:), '-', ...
       amplitudes, relativeErrors_w0{2,1}(:), '+', ...
       amplitudes, relativeErrors_w0{3,1}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,1}(:), '*', ...
       amplitudes, relativeErrors_w0{1,5}(:), '*', ...
       amplitudes, relativeErrors_w0{2,3}(:), '+', ...
       amplitudes, relativeErrors_w0{3,3}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,3}(:), '*')
legend('Regular Perturbation order 0','Regular Perturbation order 1','Regular Perturbation order 2', ...
'Regular Perturbation order 3', ...
'LP 2016 order 3', ...
'LP order 1', ...
'LP order 2', ...
'LP order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [gray; cm(3,:); cm(3,:); cm(3,:); cm(4,:);
                 cm(1,:); cm(1,:); cm(1,:)];

loglog(amplitudes, relativeErrors_w1{1,1}(:), '-', ...
       amplitudes, relativeErrors_w1{2,1}(:), '+', ...
       amplitudes, relativeErrors_w1{3,1}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,1}(:), '*', ...
       amplitudes, relativeErrors_w1{1,5}(:), '*', ...
       amplitudes, relativeErrors_w1{2,3}(:), '+', ...
       amplitudes, relativeErrors_w1{3,3}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,3}(:), '*')
legend('Regular Perturbation order 0','Regular Perturbation order 1','Regular Perturbation order 2', ...
'Regular Perturbation order 3', ...
'LP 2016 order 3', ...
'LP order 1', ...
'LP order 2', ...
'LP order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [gray; cm(3,:); cm(3,:); cm(3,:); cm(4,:);
                 cm(1,:); cm(1,:); cm(1,:)];

loglog(amplitudes, relativeErrors_w0{2,3}(:), '+', ...
       amplitudes, relativeErrors_w0{3,3}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,3}(:), '*', ...
       amplitudes, relativeErrors_w0{2,4}(:), '+', ...
       amplitudes, relativeErrors_w0{3,4}(:), 'x', ...
       amplitudes, relativeErrors_w0{4,4}(:), '*')
legend('LP order 1','LP order 2','LP order 3', ...
'LPgamma order 1', ...
'LPgamma order 2', ...
'LPgamma order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [cm(3,:); cm(3,:); cm(3,:);
                 cm(1,:); cm(1,:); cm(1,:)];

loglog(amplitudes, relativeErrors_w1{2,3}(:), '+', ...
       amplitudes, relativeErrors_w1{3,3}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,3}(:), '*', ...
       amplitudes, relativeErrors_w1{2,4}(:), '+', ...
       amplitudes, relativeErrors_w1{3,4}(:), 'x', ...
       amplitudes, relativeErrors_w1{4,4}(:), '*')
legend('LP order 1','LP order 2','LP order 3', ...
'LPgamma order 1', ...
'LPgamma order 2', ...
'LPgamma order 3', ...
'Location', 'NorthWest')
title('Bogdanov Takens normal form homoclinic convergence plot')
xlabel('$A_0$')
ylabel('$\delta(X)$')
ax = gca;
gray = [0.8 0.8 0.8];
ax.ColorOrder = [cm(3,:); cm(3,:); cm(3,:);
                 cm(1,:); cm(1,:); cm(1,:)];

writematrix([amplitudes', relativeErrors_w0{1,1}(:)], '../../data/RPorder0w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{2,1}(:)], '../../data/RPorder1w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{3,1}(:)], '../../data/RPorder2w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{4,1}(:)], '../../data/RPorder3w0.csv', 'Delimiter', ' ')
                                                   
writematrix([amplitudes', relativeErrors_w0{2,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder1w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{3,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder2w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{4,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder3w0.csv', 'Delimiter', ' ')
                                                   
writematrix([amplitudes', relativeErrors_w1{1,1}(:)], '../../data/RPorder0w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{2,1}(:)], '../../data/RPorder1w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{3,1}(:)], '../../data/RPorder2w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{4,1}(:)], '../../data/RPorder3w1.csv', 'Delimiter', ' ')
                                                   
writematrix([amplitudes', relativeErrors_w1{2,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder1w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{3,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder2w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{4,2}(:)], '../../data/RPL2PhaseConditionCorrectedorder3w1.csv', 'Delimiter', ' ')

writematrix([amplitudes(1:2:end)', relativeErrors_w0{1,5}(1:2:end)'], '../../data/LPWithoutTimeInversionorder3w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{2,3}(:)], '../../data/orbitalorder1w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{3,3}(:)], '../../data/orbitalorder2w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{4,3}(:)], '../../data/orbitalorder3w0.csv', 'Delimiter', ' ')

writematrix([amplitudes(1:2:end)', relativeErrors_w1{1,5}(1:2:end)'], '../../data/LPWithoutTimeInversionorder3w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{2,3}(:)], '../../data/orbitalorder1w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{3,3}(:)], '../../data/orbitalorder2w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{4,3}(:)], '../../data/orbitalorder3w1.csv', 'Delimiter', ' ')

writematrix([amplitudes', relativeErrors_w0{2,4}(:)], '../../data/L2PhaseConditionorder1w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{3,4}(:)], '../../data/L2PhaseConditionorder2w0.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w0{4,4}(:)], '../../data/L2PhaseConditionorder3w0.csv', 'Delimiter', ' ')

writematrix([amplitudes', relativeErrors_w1{2,4}(:)], '../../data/L2PhaseConditionorder1w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{3,4}(:)], '../../data/L2PhaseConditionorder2w1.csv', 'Delimiter', ' ')
writematrix([amplitudes', relativeErrors_w1{4,4}(:)], '../../data/L2PhaseConditionorder3w1.csv', 'Delimiter', ' ')

% example file for a single model run, with a depth dependent input
% parameter

% initialize input Scenario 1: variable friction
sc1_variable_friction = PantherAnalysis;
% change some input properties
sc1_variable_friction.setInputParameter('width_HW', 500);
sc1_variable_friction.setInputParameter('throw', 50);


% set depth varying initial stress - sinusoidal variation
shsv_default = sc1_variable_friction.getInputParameter('shsv');
sc1_variable_friction.setDepthDependentInputParameter('shsv', ones(size(sc1_variable_friction.y))*shsv_default);
%sc1_variable_friction.input_parameters.shsv.value_with_depth = ones(size(sc1_variable_friction.y))*shsv_default;
i_mid = ceil(length(sc1_variable_friction.y)/2);
% introduce a perturbation
F = 10;     % freq 
F2 = 0.9;
amp = 0.025;    % amplitude
amp2 = 0.015;
pert = amp.*sin(2*pi*F.*sc1_variable_friction.y);
sc1_variable_friction.setDepthDependentInputParameter('shsv',...
    sc1_variable_friction.getInputParameter('shsv') + pert) ;

% generate model ensemble
sc1_variable_friction.generate_ensemble();

% run panther with current input instance
sc1_variable_friction.run();

%% plot the results
h1 = figure(1); clf(h1);

subplot(1,3,1)
hold on
y_abs = sc1_variable_result.getDepth();
hp(1) = plot(sc1_variable_friction.faultResults.sne(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(2) = plot(sc1_variable_friction.faultResults.sne(:,end), y_abs);

hp(3) = plot(sc1_variable_friction.faultResults.tau(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(4) = plot(sc1_variable_friction.faultResults.tau(:,end), y_abs);

set(hp, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc1_variable_friction.getInputParameter('depth_mid') - 300, ...
    sc1_variable_friction.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hp, {'Initial normal stress', 'Final normal stress',...
    'Initial shear stress', 'Final shear stress'});

subplot(1,3,2); hold on;
scu = sc1_variable_friction.get_scu();
hs(1) = plot(scu(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(scu(:,end), y_abs);
set(hs, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc1_variable_friction.getInputParameter('depth_mid') - 300, ...
    sc1_variable_friction.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});


subplot(1,3,3); hold on;

hs(1) = plot(sc1_variable_friction.faultResults.slip(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(sc1_variable_friction.faultResults.slip(:,end), y_abs);
set(hs, 'LineWidth', 1.5);
xlabel('Slip (m)');
ylabel('Depth (m)');
ylim([sc1_variable_friction.getInputParameter('depth_mid') - 300, ...
    sc1_variable_friction.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});

%%

% initialize input Scenario 2: variable dip
sc2_variable_dip = PantherAnalysis;

% set depth varying dip
% sc2_variable_dip.input_parameters.dip.uniform_with_depth = 0;       % make dip depth-variable
y = sc2_variable_dip.y;
varying_dip = 80*ones(size(y));
varying_dip(i_mid:end) = 60;
sc2_variable_dip.setDepthDependentInputParameter('dip', varying_dip);
% turn aseismic slip off
sc2_variable_dip.aseismic_slip = 0;

% generate model ensemble
sc2_variable_dip.generate_ensemble();

% run panther with current input instance
sc2_variable_dip.run();

% plot the results Scenario 2
h2 = figure(2); clf(h2);

subplot(1,3,1)
hold on
y_abs = sc2_variable_dip.getDepth();
hp(1) = plot(sc2_variable_dip.faultResults.sne(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(2) = plot(sc2_variable_dip.faultResults.sne(:,end), y_abs);

hp(3) = plot(sc2_variable_dip.faultResults.tau(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(4) = plot(sc2_variable_dip.faultResults.tau(:,end), y_abs);

set(hp, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc2_variable_dip.getInputParameter('depth_mid') - 300, sc2_variable_dip.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hp, {'Initial normal stress', 'Final normal stress',...
    'Initial shear stress', 'Final shear stress'});

subplot(1,3,2); hold on;
scu = sc2_variable_dip.get_scu();
hs(1) = plot(scu(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(scu(:,end), y_abs);
set(hs, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc2_variable_dip.getInputParameter('depth_mid') - 300, sc2_variable_dip.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});
title('Example: Dip changes from 80 to 60 below mid reservoir');

subplot(1,3,3); hold on;

hs(1) = plot(sc2_variable_dip.faultResults.slip(:,1), y_abs, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(sc2_variable_dip.faultResults.slip(:,end), y_abs);
set(hs, 'LineWidth', 1.5);
xlabel('Slip (m)');
ylabel('Depth (m)');
ylim([sc2_variable_dip.getInputParameter('depth_mid') - 300, sc2_variable_dip.getInputParameter('depth_mid') + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});

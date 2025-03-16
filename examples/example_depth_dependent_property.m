% example file for a single model run, with a depth dependent input
% parameter

% initialize input Scenario 1: variable friction
sc1_variable_friction = PantherInput;
% change some input properties
sc1_variable_friction.input_parameters.width_HW.value = 500;
sc1_variable_friction.input_parameters.throw.value = 50;

% set depth varying initial stress - sinusoidal variation
shsv_default = sc1_variable_friction.input_parameters.shsv.value;
sc1_variable_friction.input_parameters.shsv.value_with_depth = ones(size(sc1_variable_friction.y))*shsv_default;
i_mid = ceil(length(sc1_variable_friction.y)/2);
% introduce a perturbation
F = 10;     % freq 
F2 = 0.9;
amp = 0.025;    % amplitude
amp2 = 0.015;
pert = amp.*sin(2*pi*F.*sc1_variable_friction.y);
sc1_variable_friction.input_parameters.shsv.value_with_depth = sc1_variable_friction.input_parameters.shsv.value_with_depth + pert ;
sc1_variable_friction.input_parameters.shsv.uniform_with_depth = 0;

% generate model ensemble
sc1_variable_friction.generate_ensemble();

% run panther with current input instance
sc1_variable_result = panther(sc1_variable_friction);

%% plot the results
h1 = figure(1); clf(h1);

subplot(1,3,1)
hold on
y = sc1_variable_result.y + sc1_variable_result.ensemble.depth_mid;
hp(1) = plot(sc1_variable_result.stress{1}.sne(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(2) = plot(sc1_variable_result.stress{1}.sne(:,end), y);

hp(3) = plot(sc1_variable_result.stress{1}.tau(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(4) = plot(sc1_variable_result.stress{1}.tau(:,end), y);

set(hp, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc1_variable_result.ensemble.depth_mid - 300, sc1_variable_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hp, {'Initial normal stress', 'Final normal stress',...
    'Initial shear stress', 'Final shear stress'});

subplot(1,3,2); hold on;
scu = sc1_variable_result.stress{1}.get_scu(0.6,0);
hs(1) = plot(scu(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(scu(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc1_variable_result.ensemble.depth_mid - 300, sc1_variable_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});


subplot(1,3,3); hold on;

hs(1) = plot(sc1_variable_result.slip{1}.slip(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(sc1_variable_result.slip{1}.slip(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Slip (m)');
ylabel('Depth (m)');
ylim([sc1_variable_result.ensemble.depth_mid - 300, sc1_variable_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});

%%

% initialize input Scenario 2: variable dip
sc2_variable_dip = PantherInput;

% set depth varying dip
sc2_variable_dip.input_parameters.dip.uniform_with_depth = 0;       % make dip depth-variable
y = sc2_variable_dip.y;
sc2_variable_dip.input_parameters.dip.value_with_depth = 80*ones(size(y));
i_mid = floor(length(y)/2);
sc2_variable_dip.input_parameters.dip.value_with_depth(i_mid:end) = 60;
% sc2_variable_dip.input_parameters.dip.value_with_depth = 60 + y*dip_gradient_per_m;

%sc2_variable_dip.input_parameters.width_HW.value = 0;

% turn aseismic slip off
sc2_variable_dip.aseismic_slip = 0;

% generate model ensemble
sc2_variable_dip.generate_ensemble();

% run panther with current input instance
sc2_result = panther(sc2_variable_dip);


% plot the results Scenario 2
h2 = figure(2); clf(h2);

subplot(1,3,1)
hold on
y = sc2_result.y + sc2_result.ensemble.depth_mid;
hp(1) = plot(sc2_result.stress{1}.sne(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(2) = plot(sc2_result.stress{1}.sne(:,end), y);

hp(3) = plot(sc2_result.stress{1}.tau(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(4) = plot(sc2_result.stress{1}.tau(:,end), y);

set(hp, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc2_result.ensemble.depth_mid - 300, sc2_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hp, {'Initial normal stress', 'Final normal stress',...
    'Initial shear stress', 'Final shear stress'});

subplot(1,3,2); hold on;
scu = sc2_result.stress{1}.get_scu(0.6,0);
hs(1) = plot(scu(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(scu(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([sc2_result.ensemble.depth_mid - 300, sc2_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});
title('Example: Dip changes from 90 to 60 below mid reservoir');

subplot(1,3,3); hold on;

hs(1) = plot(sc2_result.slip{1}.slip(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(sc2_result.slip{1}.slip(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Slip (m)');
ylabel('Depth (m)');
ylim([sc2_result.ensemble.depth_mid - 300, sc2_result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});

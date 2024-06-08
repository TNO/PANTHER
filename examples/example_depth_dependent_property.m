% example file for a single model run, with a depth dependent input
% parameter

% initialize input 
run_instance = PantherInput;
% change some input properties
run_instance.input_parameters.width_HW.value = 500;
run_instance.input_parameters.throw.value = 50;

% set depth varying initial stress - sinusoidal variation
shsv_default = run_instance.input_parameters.shsv.value;
run_instance.input_parameters.shsv.value_with_depth = ones(size(run_instance.y))*shsv_default;
i_mid = ceil(length(run_instance.y)/2);
% introduce a perturbation
F = 10;     % freq 
F2 = 0.9;
amp = 0.025;    % amplitude
amp2 = 0.015;
pert = amp.*sin(2*pi*F.*run_instance.y);
%pert2 = amp2.*sin(2*pi*F2*run_instance.y);
run_instance.input_parameters.shsv.value_with_depth = run_instance.input_parameters.shsv.value_with_depth + pert ;
run_instance.input_parameters.shsv.uniform_with_depth = 0;

% generate model ensemble
run_instance.generate_ensemble();
% run panther with current input instance
result = panther(run_instance);

%% plot the results
h1 = figure(1); clf(h1);

subplot(1,3,1)
hold on
y = result.y + result.ensemble.depth_mid;
hp(1) = plot(result.stress{1}.sne(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(2) = plot(result.stress{1}.sne(:,end), y);

hp(3) = plot(result.stress{1}.tau(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hp(4) = plot(result.stress{1}.tau(:,end), y);

set(hp, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([result.ensemble.depth_mid - 300, result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hp, {'Initial normal stress', 'Final normal stress',...
    'Initial shear stress', 'Final shear stress'});

subplot(1,3,2); hold on;
scu = result.stress{1}.get_scu(0.6,0);
hs(1) = plot(scu(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(scu(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Stress (MPa)');
ylabel('Depth (m)');
ylim([result.ensemble.depth_mid - 300, result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});


subplot(1,3,3); hold on;

hs(1) = plot(result.slip{1}.slip(:,1), y, 'LineStyle','--','Color',[0.5,0.5,0.5]);
hs(2) = plot(result.slip{1}.slip(:,end), y);
set(hs, 'LineWidth', 1.5);
xlabel('Slip (m)');
ylabel('Depth (m)');
ylim([result.ensemble.depth_mid - 300, result.ensemble.depth_mid + 300]);
set(gca,'Box',1);
legend(hs, {'Initial', 'Final'});

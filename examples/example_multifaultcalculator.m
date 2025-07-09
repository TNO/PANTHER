% example for using Multifaultcalculator
% Multifaultcalculator handles multiple 2D fault cross-sections (pillars)

% read example fault input
fault_segments = readtable(fullfile(fileparts(matlab.desktop.editor.getActiveFilename), 'example_files\example_fault_input_for_multifaultcalculator'));

% create MultiFaultCalculator object
fault = FaultSurfaceCalculator(height(fault_segments));

% save the metadata of interest, which will be added to the fault info in
% the fault instance of MultiFaultCalculator
metadata_to_add = fault_segments(:,2:end);
fault = fault.add_pillar_info_as_table(metadata_to_add);

% specify an example net-to-gross grid and add to pillar metadata
x_vector = linspace(min(fault.pillar_info.X), max(fault.pillar_info.X), 12);
y_vector = linspace(min(fault.pillar_info.Y), max(fault.pillar_info.Y), 12);
[X_grid, Y_grid] = meshgrid(x_vector, y_vector);
ntg_grid = ones(size(X_grid))*0.8;
ntg_grid = ntg_grid + repmat(linspace(0,10,size(X_grid,2))*0.02,size(X_grid,1), 1);
fault = fault.add_info_from_closest_point(X_grid, Y_grid, ntg_grid, 'NTG');
friction_from_ntg = 0.6 - 0.2*(1-fault.pillar_info.NTG);  % arbitrary example relation NTG and friction
clear x_vector y_vector X_grid Y_grid ntg_grid

% generate example of a depth dependent property
y = fault.pillars{1}.y;
depth_dependent_shsv = cell(size(fault.pillars));
depth_dependent_shsv(:) = {ones(size(y))*0.8};
% set the reservoir interval shsv to a different value
for i = 1 : length(fault.pillars)
    fault.pillars{i}.generate_ensemble();
    top_HW_i = fault.pillars{i}.ensemble_members{1}.top_HW_i(y);
    top_FW_i = fault.pillars{i}.ensemble_members{1}.top_FW_i(y);
    base_HW_i = fault.pillars{i}.ensemble_members{1}.base_HW_i(y);
    base_FW_i = fault.pillars{i}.ensemble_members{1}.base_FW_i(y);
    top_reservoir_interval = min(top_FW_i, top_HW_i);
    base_reservoir_interval = max(base_FW_i, base_HW_i);
    depth_dependent_shsv{i}(top_reservoir_interval:base_reservoir_interval) = 0.75;
end

% change the input for the different fault pillars
fault.set_input_parameter('dip', fault_segments.dip);
fault.set_input_parameter('dip_azi', fault_segments.azimuth);
fault.set_input_parameter('f_s', friction_from_ntg);
fault.set_depth_dependent_input_parameter('shsv', depth_dependent_shsv);

% change the run settings (pillar settings like diffusion_P, save_stress,
% aseismic_slip, nucleation_criterion, etc. )
fault.set_run_setting('diffusion_P', 1);

% change MultiFaultCalculator settings
fault.parallel = 0;
fault.suppress_pillar_run_status_output  = 0;

% change pressure loading on the last 5 pillars
load_tables = cell(height(fault_segments), 1);
load_tables(:) = {initialize_load_table()};     % initialize default load table
for j = 6 : height(fault_segments)
    % set hanging wall depletion to zero
    load_tables{j}.P_factor_HW = zeros(height(load_tables{j}.P_steps),1);
end
fault = fault.set_load_tables(load_tables);

% run the simulation over all pillars of the fault
fault = fault.run();
fault = fault.reduce_output([1:5:35]);

% plot some summary results 
h1 = figure(1); clf(h1);
subplot(3,1,1)
hold on
plot(fault.pillar_info.Y,fault.result_summary.reactivation_dP);
plot(fault.pillar_info.Y,fault.result_summary.nucleation_dP);
legend({'Reactivation pressure change (MPa)','Nucleation pressure change (MPa)'});
xlabel('Distance along y');
ylabel('\DeltaP (MPa)');
subplot(3,1,2)
hold on
plot(fault.pillar_info.Y,fault.pillar_info.dip);
xlabel('Distance along Y');
ylabel('Dip (deg)');
subplot(3,1,3)
hold on
plot(fault.pillar_info.Y,fault.pillar_info.azimuth);
legend({ 'Dip azimuth'});
xlabel('Distance along Y');
ylabel('Dip azimuth (deg)');

%% example plot of an input parameter along the fault surface
[grid_dip, along_fault_length_grid, depth_grid, ~] = fault.get_fault_grid_for_input_parameter('dip');
[grid_shsv] = fault.get_fault_grid_for_input_parameter('shsv');
[tau_on_fault] = fault.get_fault_grid_for_output('tau', 7);
[sne_on_fault] = fault.get_fault_grid_for_output('sne', 7);

tops = fault.get_top_base_reservoir();

h2 = figure(2); clf(h2);

subplot(3,1,1)
hold on
hs = surf(along_fault_length_grid, depth_grid, grid_dip);
view([0,90]);
set(hs, 'EdgeColor','none');
xlabel('Along fault length (m)');
ylabel('Depth (m)');
cb = colorbar();
ylabel(cb, 'Dip');
xlim([min(fault.L_strike), max(fault.L_strike)]);
[ymin, ymax]  = fault.get_min_max_depth;
ylim([ymin, ymax]);
plot3(fault.L_strike, tops.top_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.base_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.top_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 
plot3(fault.L_strike, tops.base_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 

subplot(3,1,2)
hold on
hs = surf(along_fault_length_grid, depth_grid, grid_shsv);
view([0,90]);
set(hs, 'EdgeColor','none');
xlabel('Along fault length (m)');
ylabel('Depth (m)');
cb = colorbar();
ylabel(cb, '\sigma_h / \sigma_v');
xlim([min(fault.L_strike), max(fault.L_strike)]);
[ymin, ymax]  = fault.get_min_max_depth;
ylim([ymin, ymax]);
plot3(fault.L_strike, tops.top_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.base_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.top_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 
plot3(fault.L_strike, tops.base_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 

subplot(3,1,3)
hold on
hs = surf(along_fault_length_grid, depth_grid, tau_on_fault./sne_on_fault);
view([0,90]);
set(hs, 'EdgeColor','none');
xlabel('Along fault length (m)');
ylabel('Depth (m)');
cb = colorbar();
ylabel(cb, '\tau/\sigma_n''');
xlim([min(fault.L_strike), max(fault.L_strike)]);
[ymin, ymax]  = fault.get_min_max_depth;
ylim([ymin, ymax]);
plot3(fault.L_strike, tops.top_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.base_FW, ones(size(fault.L_strike))*1e8, 'LineStyle','-','Color','k'); 
plot3(fault.L_strike, tops.top_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 
plot3(fault.L_strike, tops.base_HW, ones(size(fault.L_strike))*1e8, 'LineStyle','--','Color','k'); 



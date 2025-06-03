% example for using Multifaultcalculator
% Multifaultcalculator handles multiple 2D fault cross-sections (pillars)

% read example fault input
fault_segments = readtable(fullfile(fileparts(matlab.desktop.editor.getActiveFilename), 'example_files\example_fault_input_for_multifaultcalculator'));

% create MultiFaultCalculator object
fault = MultiFaultCalculator(height(fault_segments));

% save the metadata of interest, which will be added to the fault info in
% the fault instance of MultiFaultCalculator
metadata_to_add = fault_segments(:,2:end);
fault = fault.add_pillar_info_as_table(metadata_to_add);

% specify an example net-to-gross grid and add to pillar metadata
x_vector = linspace(min(fault.pillar_info.X), max(fault.pillar_info.Y), 12);
y_vector = linspace(min(fault.pillar_info.Y), max(fault.pillar_info.Y), 12);
[X_grid, Y_grid] = meshgrid(x_vector, y_vector);
ntg_grid = ones(size(X_grid))*0.8;
ntg_grid = ntg_grid + repmat(linspace(0,10,size(X_grid,2))*0.02,size(X_grid,1), 1);
fault = fault.add_info_from_closest_point(X_grid, Y_grid, ntg_grid, 'NTG');
friction_from_ntg = 0.6 - 0.2*(1-fault.pillar_info.NTG);  % arbitrary example relation NTG and friction
clear x_vector y_vector X_grid Y_grid ntg_grid

% change the input for the different fault pillars
fault.set_input_parameter('dip', fault_segments.dip);
fault.set_input_parameter('dip_azi', fault_segments.azimuth);
fault.set_input_parameter('f_s', friction_from_ntg);

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

% plot the fault surface
[along_fault_length_grid, depth_grid, values, ~] = fault.get_fault_grid_for_input_parameter('dip');







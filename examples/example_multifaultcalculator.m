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

% change the input for the different fault pillars
fault.set_input_parameter('dip', fault_segments.dip);
fault.set_input_parameter('dip_azi', fault_segments.azimuth);

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







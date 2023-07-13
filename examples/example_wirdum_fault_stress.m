% example function for 2.5D case
% computes stresses and stress changes for each pillar of the Wirdum fault
% assumes the 2D plane-strain approximation holds at each fault pillar
% Wirdum fault file is found in panther-temp/examples/example_files

% read Wirdum fault properties
filePath = matlab.desktop.editor.getActiveFilename;
example_folder = [fileparts(filePath),'/'];
cd(example_folder);
wirdum = readtable('example_files/Wirdum_fault_reservoir_geometry_RD.csv');
wirdum = renamevars(wirdum, ["dip_azimuth","dip_angle","cdepth"], ["dip_azi","dip","depth_mid"]);

% initialize input instance
analysis = PantherInput;
% set some input parameters
analysis.input_parameters.sH_dir.value = 140;
analysis.input_parameters.shsv.value = 0.75;
analysis.diffusion_P = 0;
% analysis.aseismic_slip = 0; 
% load table values to generate ensemble
analysis.generate_ensemble_from_table(wirdum);

% running analysis (starting up parallel can take a while if it is not yet
% running)
analysis_result = panther(analysis);

% time/load index of nucleation
nucleation = min(analysis_result.summary.nucleation_index);

% interpolate stresses at the nucleation time index, for each pillar
pillar_stress = cell(height(analysis_result.summary), 1);
if ~isnan(nucleation)
    for i = 1 : length(analysis_result.stress)
        pillar_stress{i} = table();
        pillar_stress{i}.y = analysis_result.y + analysis.ensemble{i}.depth_mid;
        pillar_stress{i}.x = ones(size(pillar_stress{i}.y))*wirdum.x_coor(i);
        [pillar_stress{i}.sne, pillar_stress{i}.tau] = analysis_result.stress{i}.get_failure_stress(nucleation);
        % writetable(pillar_stress{i}, ['pillar', num2str(i)]);
    end
end

% get initial stress and stress changes with respect to initial time step
[sne0, tau0] = analysis_result.get_initial_stress();
[dsne, dtau] = analysis_result.get_stress_changes();
%%
pillarnr=50;
subplot(1,6,1);
plot(sne0{pillarnr,1},analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))

subplot(1,6,2);
plot(tau0{pillarnr,1},analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))

subplot(1,6,3);
plot(dsne{pillarnr,1}(:,35),analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))

subplot(1,6,4);
plot(dtau{pillarnr,1}(:,35),analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))

subplot(1,6,5);
plot(sne0{pillarnr,1}+dtau{pillarnr,1}(:,35),analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))

subplot(1,6,6);
plot(tau0{pillarnr,1}+dtau{pillarnr,1}(:,35),analysis_result.y+analysis_result.ensemble.depth_mid(pillarnr))




%%
% results=[sne0, tau0, dsne, dtau];
% write results to .csv file
write_2Dstress_to_csv(wirdum,analysis_result,'./output/','_diff')
% example function for a 2.5D field application, case Wirdum fault
% computes stresses and stress changes for each pillar of the Wirdum fault
% assumes the 2D plane-strain approximation holds at each fault pillar
% Wirdum fault file is found in panther-temp/examples/example_files
openProject('./../Panther.prj');

%% ------------READ FIELD DATA--------------------------
% Open the log file for writing
faultname='gocad';
% read Wirdum fault properties
runcase='gasfilled_min_dc_scaled10';

rootout=['/quanta1/home/vheiden/work/PANTHER/output_files/' faultname '/' runcase '/'];
outdir=strcat(rootout,'/');
mkdir(outdir);

% Get the current date and time
current_time = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
% Generate the log file name with the timestamp
log_file_name = sprintf('PANTHER_%s.out', datestr(current_time, 'yyyymmdd_HHMMSS'));
log_file = fopen(fullfile(outdir, log_file_name), 'w');
% Check if the file was opened successfully
if log_file == -1
    error('Failed to open log file');
end


% Retrieve the current Git branch name
[status, branch_name] = system('git rev-parse --abbrev-ref HEAD');
if status == 0
    git_branch = strtrim(branch_name);
else
    git_branch = 'Failed to retrieve Git branch name';
end

% Write the Git branch name and timestamp to the log file
fprintf(log_file, 'Git branch: %s\n', git_branch);
fprintf(log_file, 'Date and time: %s\n\n', datestr(current_time, 'yyyy-MM-dd HH:mm:ss'));



wirdum = readtable(['./../input_files/' faultname filesep 'initial_stresses_wirdum_' runcase '.csv']);
depth_properties = wirdum(:,8:end);
wirdum(:,8:end) = [];

% Write parameters to the log file
fprintf(log_file, 'fault geometry: %s\n', faultname);
fprintf(log_file, 'runecase: %s\n', runcase);


% read pressures around the Wirdum fault
pressure_file='XY_PRF_Avg_OS1_CY.csv';
pres = readtable(['./examples/example_files' filesep pressure_file]);
for i = 1 : height(wirdum)
    % find nearest neighbour    
    [~, i_min_x] = min(abs(wirdum.x_coor(i) - pres.X));
    i_min_x_all = find(pres.X == pres.X(i_min_x));
    [~, i_min_y] = min(abs(wirdum.y_coor(i) - pres.Y(i_min_x_all)));
    i_min = i_min_x_all(i_min_y);
    wirdum.pressure_years(i) = {linspace(1957, 2052, 96)'};
    wirdum.pressure_values(i) = {(table2array(pres(i_min,3:end))/10)'};
end
fprintf(log_file, 'reading pressures from: %s\n', pressure_file);

% Wirdum seismic swarm
eq_location_file='foreshock_locations_Wirdum';
wirdum_eq = readtable(['./examples/example_files' filesep eq_location_file]);
wirdum_eq.date = cellfun(@(x) x(1:10), wirdum_eq.origt,'UniformOutput',false);
wirdum_eq.years_production = datenum(wirdum_eq.date - datetime('1957-01-01'))/365.25;
wirdum_eq.years_from_first_event = wirdum_eq.years_production - wirdum_eq.years_production(1);
for i = 1 : height(wirdum_eq)
    % find nearest neighbour  
    [~, i_min_x] = min(abs(wirdum_eq.rdx(i) - pres.X));
    i_min_x_all = find(pres.X == pres.X(i_min_x));
    [~, i_min_y] = min(abs(wirdum_eq.rdy(i) - pres.Y(i_min_x_all)));
    i_min = i_min_x_all(i_min_y);
    pressure_years = linspace(1957, 2052, 96);
    pressure_values = (table2array(pres(i_min,3:end))/10);
    wirdum_eq.pressure_at_eq_time(i) = interp1(pressure_years, pressure_values, wirdum_eq.years_production(i) + 1957);
    wirdum_eq.pressure_change_at_eq_time(i) = wirdum_eq.pressure_at_eq_time(i) - pressure_values(1);

    dist_to_fault = ((wirdum.x_coor - wirdum_eq.rdx(i)).^2 + (wirdum.y_coor - wirdum_eq.rdy(i)).^2).^0.5;
    [wirdum_eq.dist_to_fault(i), i_min] = min(dist_to_fault);
    wirdum_eq.fault_rdx(i) = wirdum.x_coor(i_min);
    wirdum_eq.fault_rdy(i) = wirdum.y_coor(i_min);
end

fprintf(log_file, 'reading seismic swarm from: %s\n', eq_location_file);

%% ----------------------PREPARE AND RUN Model SIMULATION-------------

wirdum_variables = {'depth_mid','dip','dip_azi','throw','thick'};
depth_variable = {'f_s','f_d','d_c','sv_grad','shsv'};
Zechstein_defaults = [1, 0.4, 0.06, 22, 1];
for i = 1 : height(wirdum)
    pillar{i} = PantherInput;
    % set some input parameters
    pillar{i}.diffusion_P = 1;
    pillar{i}.aseismic_slip = 1; 
    pillar{i}.input_parameters.sH_dir.value = 140;
    pillar{i}.input_parameters.sHsh.value = 1.07;        % from Mmax report, and Eijs (2015)
    pillar{i}.input_parameters.p_grad.value = 10.2;      % [MPa/km] from initial stress study
    pillar{i}.input_parameters.p_grad_res.value = 0.2;   % [MPa/km] pressure gradient in reservoir
    pillar{i}.input_parameters.p_over.value = 3;         % [MPa] pressure gradient in reservoir
    pillar{i}.input_parameters.biot.value = 1;
    pillar{i}.input_parameters.poisson.value = 0.2;
    pillar{i}.nucleation_criterion = 'Day3D';
    % pillar{i}.nucleation_criterion = 'UR2D';

    for j = 1 : length(wirdum_variables)
        var = wirdum_variables{j};
        pillar{i}.input_parameters.(var).value = wirdum.(var)(i);
    end
    for j = 1 : length(depth_variable)
        var = depth_variable{j};
        default_value = Zechstein_defaults(j);  % assign value for ZE as default
        % initialize depth array for depth-dependent variables
        pillar{i}.input_parameters.(var).value_with_depth = ones(size(pillar{i}.y))*default_value;
        pillar{i}.input_parameters.(var).uniform_with_depth = 0;
        y_abs = pillar{i}.y + pillar{i}.input_parameters.depth_mid.value;


        % fill depth array, depending on horizon depths
        for k = 1 : 7
            below_horizon =  y_abs < depth_properties.(['h',num2str(k)])(i);
            pillar{i}.input_parameters.(var).value_with_depth(below_horizon) = depth_properties.([var,num2str(k)])(i);
        end
    end
    pillar{i}.generate_ensemble;
    pillar_result{i} = panther(pillar{i});
    if i == 1
        input_table = pillar{i}.ensemble_to_table();
        analysis_result = pillar_result{i};
    else
        input_table = vertcat(input_table, pillar{i}.ensemble_to_table);
        analysis_result.pressure{1,i} = pillar_result{i}.pressure{1};
        analysis_result.temperature{1,i} = pillar_result{i}.temperature{1};
        analysis_result.stress{1,i} = pillar_result{i}.stress{1};
        analysis_result.slip{1,i} = pillar_result{i}.slip{1};
        analysis_result.summary = vertcat(analysis_result.summary, pillar_result{i}.summary);
                    analysis_result.ensemble = vertcat(analysis_result.ensemble, pillar_result{i}.ensemble);
                end
            end
i=1
fprintf(log_file, '  sH_dir: %.2f\n', pillar{i}.input_parameters.sH_dir.value);
fprintf(log_file, '  sHsh: %.2f\n', pillar{i}.input_parameters.sHsh.value);
fprintf(log_file, '  p_grad: %.2f MPa/km\n', pillar{i}.input_parameters.p_grad.value);
fprintf(log_file, '  p_grad_res: %.2f MPa/km\n', pillar{i}.input_parameters.p_grad_res.value);
fprintf(log_file, '  p_over: %.2f MPa\n', pillar{i}.input_parameters.p_over.value);
fprintf(log_file, '  biot: %.2f\n', pillar{i}.input_parameters.biot.value);
fprintf(log_file, '  poisson: %.2f\n', pillar{i}.input_parameters.poisson.value);
fprintf(log_file, '  nucleation_criterion: %s\n', pillar{i}.nucleation_criterion);

% Loop through each depth-dependent variable
fprintf(log_file, 'Depth-dependent variables:\n');
for j = 1:numel(depth_variable)
    fprintf(log_file, '  %s\n', depth_variable{j});
end

% time/load index of nucleation
nucleation = min(analysis_result.summary.nucleation_index);
reactivation = min(analysis_result.summary.reactivation_index);

fprintf(log_file, 'nucleation at : %s MPa\n', mat2str(nucleation));
fprintf(log_file, 'reactivation at : %s MPa\n', mat2str(reactivation));

nucleation
reactivation

%%
% interpolate stresses at the nucleation time index, for each pillar
pillar_stress = cell(height(analysis_result.summary), 1);
if ~isnan(nucleation)
    for i = 1 : length(analysis_result.stress)
        pillar_stress{i} = table();
        pillar_stress{i}.y = analysis_result.y + analysis_result.ensemble.depth_mid(i);
        pillar_stress{i}.x = ones(size(pillar_stress{i}.y))*wirdum.x_coor(i);
        [pillar_stress{i}.sne, pillar_stress{i}.tau] = analysis_result.stress{i}.get_failure_stress(nucleation);
    end
end
%%
% get initial stress and stress changes with respect to initial time step
[sne0, tau0] = analysis_result.get_initial_stress();
[dsne, dtau] = analysis_result.get_stress_changes();
[p_init] = analysis_result.get_initial_pressure();

%% plotting results
nucstr = sprintf('%.2f', nucleation);
plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, [ outdir runcase '_' nucstr 'MPa' '_'], 'x')
plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, [ outdir runcase '_' nucstr 'MPa' '_'], 'y')
% plot_mapView_fault_events(wirdum, wirdum_eq, './')

%% write to results to .csv
write_2Dstress_to_csv(wirdum,analysis_result,round(nucleation),[ outdir runcase '_' nucstr 'MPa' '_'],'')

% Close the log file
fclose(log_file);
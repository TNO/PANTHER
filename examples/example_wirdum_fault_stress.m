% example function for a 2.5D field application, case Wirdum fault
% computes stresses and stress changes for each pillar of the Wirdum fault
% assumes the 2D plane-strain approximation holds at each fault pillar
% Wirdum fault file is found in panther-temp/examples/example_files
openProject('./../Panther.prj');

%% ------------READ FIELD DATA--------------------------
% read Wirdum fault properties
filePath = matlab.desktop.editor.getActiveFilename;
example_folder = [fileparts(filePath),filesep];
cd(example_folder);

runcase='gasfilled_dc_scaled10';
scaled=10;
wirdum = readtable(['example_files' filesep 'initial_stresses_wirdum_' runcase '.csv']);
% wirdum = readtable(['example_files' filesep 'initial_stresses_wirdum.csv']);
depth_properties = wirdum(:,8:end);
wirdum(:,8:end) = [];

% read pressures around the Wirdum fault
pres = readtable(['example_files' filesep 'XY_PRF_Avg_OS1_CY.csv']);
for i = 1 : height(wirdum)
    % find nearest neighbour    
    [~, i_min_x] = min(abs(wirdum.x_coor(i) - pres.X));
    i_min_x_all = find(pres.X == pres.X(i_min_x));
    [~, i_min_y] = min(abs(wirdum.y_coor(i) - pres.Y(i_min_x_all)));
    i_min = i_min_x_all(i_min_y);
    wirdum.pressure_years(i) = {linspace(1957, 2052, 96)'};
    wirdum.pressure_values(i) = {(table2array(pres(i_min,3:end))/10)'};
end

% Wirdum seismic swarm
wirdum_eq = readtable(['example_files' filesep 'all_onfault_earthquakes.xlsx']);
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

%% ----------------------PREPARE AND RUN Model SIMULATION-------------

wirdum_variables = {'depth_mid','dip','dip_azi','throw','thick'};
depth_variable = {'f_s','f_d','d_c','sv_grad','shsv'};
% Zechstein_defaults = [0.8, 0.4, 0.001, 22, 0.99];
Zechstein_defaults = [1, 0.4, 0.006/scaled, 22, 1];

% Initialize the muds table with NaN values
fss = nan(501,height(wirdum));
fds = nan(501,height(wirdum));
dcs = nan(501,height(wirdum));
sv_grads = nan(501,height(wirdum));
shsvs = nan(501,height(wirdum));

for i = 1 : height(wirdum)
    pillar{i} = PantherInput;
    % set some input parameters
    pillar{i}.diffusion_P = 1;
    pillar{i}.aseismic_slip = 0; 
    pillar{i}.input_parameters.sH_dir.value = 140;
    % pillar{i}.input_parameters.shsv.value = 0.75;        % from initial stress study
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

    pill=pillar{i}.generate_ensemble;
    % Extract the depth dependent values for the current pillar
    depth_dependent_values = pillar{i}.generate_ensemble.ensemble{1};
    
    % Store the depth dependent values
    fss(:, i) = depth_dependent_values.f_s;
    fds(:, i) = depth_dependent_values.f_d;
    dcs(:, i) = depth_dependent_values.d_c;
    sv_grads(:, i) = depth_dependent_values.sv_grad;
    shsvs(:, i) = depth_dependent_values.shsv;

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



%%
% time/load index of nucleation
[nucleation, nuc_index] = min(analysis_result.summary.nucleation_index);
[reactivation, reac_index] = min(analysis_result.summary.reactivation_index);

%%
% interpolate stresses at the nucleation time index, for each pillar
pillar_stress = cell(height(analysis_result.summary), 1);
if ~isnan(nucleation)
    for i = 1 : length(analysis_result.stress)
        pillar_stress{i} = table();
        pillar_stress{i}.y = analysis_result.y + analysis_result.ensemble.depth_mid(i);
        pillar_stress{i}.x = ones(size(pillar_stress{i}.y))*wirdum.x_coor(i);
        [pillar_stress{i}.sne, pillar_stress{i}.tau] = analysis_result.stress{i}.get_failure_stress(nucleation);
        % writetable(pillar_stress{i}, ['pillar', num2str(i)]);
    end
end
%%
% get initial stress and stress changes with respect to initial time step
[sne0, tau0] = analysis_result.get_initial_stress();
[dsne, dtau] = analysis_result.get_stress_changes();
[p_init] = analysis_result.get_initial_pressure();

%% plotting results
rootout=['/Users/4146751/surfdrive2/PhD/DeepNL_Phys_Mmax/WP1/Code/PANTHER_output/output/basecase_' runcase '/'];
faultname='gocad';

outdir=strcat(rootout,'/',faultname,'/');
mkdir(outdir);
nucstr = sprintf('%.2f', nucleation);
plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, [ outdir runcase '_' nucstr 'MPa' '_'], 'x',{fss,fds,dcs,sv_grads,shsvs},round(nucleation),nuc_index,round(reactivation),reac_index)
% plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, [ outdir runcase '_' nucstr 'MPa' '_'], 'y',{fss,fds,dcs,sv_grads,shsvs},round(nucleation),nuc_index,round(reactivation),reac_index)
% plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, [ outdir runcase '_' nucstr 'MPa' '_'], 'dist',{fss,fds,dcs,sv_grads,shsvs},round(nucleation),nuc_index,round(reactivation),reac_index)

%%
% Add the directory containing the function to MATLAB's search path
% data_folder='/Users/4146751/surfdrive2/PhD/DeepNL_Phys_Mmax/WP1/Code/PANTHER_output/output/basecase_gasfilled_dc_scaled10/';
% interpolate_and_write_to_netcdf(wirdum,analysis_result,{fss,fds,dcs,sv_grads,shsvs},round(nucleation),nuc_index,round(reactivation),reac_index,strcat(outdir,faultname,'_'),'');


%%
% write results to .csv file
% 
% mkdir(outdir);
% write_2Dstress_to_csv(wirdum,analysis_result,round(nucleation),strcat(outdir,faultname,'_'),'')
% examples script to use GroningenFaultCalculator
% TODO make into a function which uses a config file

%% prepare input
proj = matlab.project.rootProject;
proj_folder = convertStringsToChars(proj.RootFolder);
% out_folder = [proj_folder,'/results/test_single_fault_segment/'];

file_name.faults = [proj_folder,'/examples/Mmax/Mmax_files/Bourne_2017_FaultModel_Geometries_DuplicatesRemoved.xlsx'];
file_name.poro = [proj_folder,'/examples/Mmax/Mmax_files/Petrel_Slochteren_POR_mean_no_header.txt'];
file_name.h_ROCLT = [proj_folder,'/examples/Mmax/Mmax_files/Petrel_thickness_ROCLT.csv'];

bourne_faults = readtable(file_name.faults);
poro = readtable(file_name.poro);
poro = renamevars(poro, {'Var1','Var2','Var3'},{'X','Y','poro'});
h_ROCLT = readtable(file_name.h_ROCLT);

fault_names = unique(bourne_faults.Fault);

%% prepare the simulation. runs per Groningen fault (fault{j})
% run settings
diffusion = 0;      % pore pressure diffusion

% number of faults to be used in the calculation 
n_faults = 2;        % set to length(fault_names) for all faults


for j = 827 : length(fault_names) %n_faults   
    
    i_fault = strcmp(bourne_faults.Fault, fault_names{j,1});

    % initialize multi fault object
    fault{j,1} = MultiFaultCalculator(sum(i_fault));
    
    % populate with input from the fault and geological model
    fault{j,1} = fault{j,1}.add_pillar_info_as_table(bourne_faults(i_fault,["Fault","Easting","Northing"]));
    fault{j,1} = fault{j,1}.set_input_parameter('dip', bourne_faults.Dip(i_fault));
    fault{j,1} = fault{j,1}.set_input_parameter('dip_azi', bourne_faults.DipAzimuth(i_fault));
    
    % specify and add thickness Slochteren Formation (ROSL). Ten Boer not depleting
    thickness = (bourne_faults.ThicknessLeft2(i_fault) + bourne_faults.ThicknessRight2(i_fault))/2;   % thickness of entire RO Group (ROSL + ROCLT)
    fault{j,1} = fault{j,1}.add_info_from_closest_point(h_ROCLT.X, h_ROCLT.Y, h_ROCLT.Thickness, 'Easting','Northing','h_ROCLT');
    fault{j,1} = fault{j,1}.add_pillar_info_as_table(table(thickness - fault{j,1}.pillar_info.h_ROCLT,'VariableNames',{'h_ROSL'}));
    fault{j,1} = fault{j,1}.set_input_parameter('thick', fault{j,1}.pillar_info.h_ROSL);
    
    % add porosity values
    fault{j,1} = fault{j,1}.add_info_from_closest_point(poro.X, poro.Y, poro.poro, 'Easting','Northing','poro');
    
    % specify the run settings (can be same or different per pillar
    fault{j,1} = fault{j,1}.set_run_setting('diffusion_P',0);
    fault{j,1} = fault{j,1}.set_run_setting('load_case',repmat({'P'},1));
    fault{j,1} = fault{j,1}.set_run_setting('aseismic_slip', 0);
    fault{j,1} = fault{j,1}.set_run_setting('nucleation_criterion', {'UR2D'});
    
    
    fault{j,1}.parallel = 1;

    % execute run
    fault{j,1} = fault{j,1}.run();
    
    nucleation_load_step = fault{j,1}.get_minimum_nucleation_load_step();
    % overwrite the pillar-specific nucleation stress with the stress at the 
    % minimum nucleation load step over the whole fault
    if ~isnan(nucleation_load_step) 
        fault{j,1} = fault{j,1}.overwrite_nucleation_stress(nucleation_load_step);
    end

    % reduce output of stress, slip, pressure, temperature arrays
    % these were not reduced earlier, because we needed to determine the
    % minimum nucleation pressure on the whole fault, and compute the
    % stresses at that nucleation pressure. for that, we need all time
    % steps. 
    fault{j,1} = fault{j,1}.reduce_output([1, length(fault{j,1}.pillars{1}.load_case)]);

end


% @Huihui stresses for each fault pillar can be found in the pillar
% results, stress object. E.g fault 2, 4th pillar, stress realization
% fault{2}.pillar_results{4}.stress{1}
% this contains the reduced sne and tau output (length(y), n_time_steps),
% as well as the sne and tau at reactivation and nucleation
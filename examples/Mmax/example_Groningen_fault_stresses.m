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


%% prepare the simulation
% run settings
diffusion = 0;      % pore pressure diffusion

% number of fault pillars to be used in the calculation 
n_faults = 1000;        % set to height of bourne_faults for all faults

% initialize multi fault object
groningen = GroningenFaultCalculator(n_faults);

% populate with input from the fault and geological model
groningen = groningen.add_pillar_info_as_table(bourne_faults(1:n_faults,["Fault","Easting","Northing"]));
groningen = groningen.set_input_parameter('dip', bourne_faults.Dip);
groningen = groningen.set_input_parameter('dip_azi', bourne_faults.DipAzimuth);

% specify and add thickness Slochteren Formation (ROSL). Ten Boer not depleting
thickness = (bourne_faults.ThicknessLeft2(1:n_faults) + bourne_faults.ThicknessRight2(1:n_faults))/2;   % thickness of entire RO Group (ROSL + ROCLT)
groningen = groningen.add_info_from_closest_point(h_ROCLT.X, h_ROCLT.Y, h_ROCLT.Thickness, 'Easting','Northing','h_ROCLT');
groningen = groningen.add_pillar_info_as_table(table(thickness - h_ROCLT.Thickness(1:n_faults),'VariableNames',{'h_ROSL'}));
groningen = groningen.set_input_parameter('thick', groningen.pillar_info.h_ROSL);

% add porosity values
groningen = groningen.add_info_from_closest_point(poro.X, poro.Y, poro.poro, 'Easting','Northing','poro');

% specify the run settings (can be same or different per pillar
groningen = groningen.set_run_setting('diffusion_P',0);
groningen = groningen.set_run_setting('load_case',repmat({'P'},1));
groningen = groningen.set_run_setting('aseismic_slip', 0);
groningen = groningen.set_run_setting('nucleation_criterion', {'UR2D'});

% execute run
groningen = groningen.run();

run_summary = groningen.get_results_summary();  % get summary of all run results

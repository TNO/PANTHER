% example function for a 2.5D field application, case Wirdum fault
% computes stresses and stress changes for each pillar of the Wirdum fault
% assumesgit p the 2D plane-strain approximation holds at each fault pillar
% Wirdum fault file is found in panther-temp/examples/example_files


%% ------------READ FIELD DATA--------------------------
% read Wirdum fault properties
filePath = matlab.desktop.editor.getActiveFilename;
example_folder = [fileparts(filePath),'\'];
cd(example_folder);
wirdum = readtable('example_files\Wirdum_fault_reservoir_geometry_RD_realistic.csv');
wirdum = renamevars(wirdum, ["dip_azimuth","dip_angle","cdepth"], ["dip_azi","dip","depth_mid"]);
wirdum.dip_azi = wirdum.dip_azi - 90;

% read pressures around the Wirdum fault
pres = readtable('example_files\XY_PRF_Avg_OS1_CY.csv');
for i = 1 : height(wirdum)
    % find nearest neighbour    
    [~, i_min_x] = min(abs(wirdum.x_coor(i) - pres.X));
    i_min_x_all = find(pres.X == pres.X(i_min_x));
    [~, i_min_y] = min(abs(wirdum.y_coor(i) - pres.Y(i_min_x_all)));
    i_min = i_min_x_all(i_min_y);
    wirdum.pressure_years(i) = {linspace(1957, 2052, 96)'};
    wirdum.pressure_values(i) = {(table2array(pres(i_min,3:end))/10)'};
end
h1 = figure(1); clf(h1);hold on
for i = 1 : height(wirdum)
    plot(wirdum.pressure_years{i}, wirdum.pressure_values{i});
end
xlabel('Year'); ylabel('Reservoir pressure at pillar (MPa)');

% Wirdum seismic swarm
wirdum_eq = readtable("example_files\foreshock_locations_Wirdum.xlsx");
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

wirdum_eq.fault_rdx


%% ----------------------PREPARE AND RUN Model SIMULATION-------------

% initialize input instance
analysis = PantherInput;

% set some run settings
analysis.diffusion_P = 1;
% analysis.aseismic_slip = 0; 

% set some input parameters
analysis.input_parameters.sH_dir.value = 140;
analysis.input_parameters.shsv.value = 0.75;
analysis.input_parameters.sHsh.value = 1.1;

% load table values to generate ensemble
% note: this assumed similar pressure drop along the fault
analysis.generate_ensemble_from_table(wirdum);
input_table = analysis.ensemble_to_table();

% running analysis (starting up parpool can take 30s when not yet running)
analysis_result = panther(analysis);

% here add a routing to subtract stress

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

fig_path = '\\tsn.tno.nl\data\sv\sv-053185\Kluis\Research\DeepNL\PhysMax\Widrum_nucleation_modeling\';
f_name = ['Wirdum_shsh',num2str(analysis.input_parameters.sHsh.value*10,'%.0f'),...
    '_diff',num2str(analysis.diffusion_P)];
save([fig_path, f_name],'wirdum','wirdum_eq','analysis','analysis_result','input_table');

%% plot dip, strike, throw, reactivation and nucleation pressure along fault with events
h2 = figure(2); clf(h2);
set(gcf, 'Units', 'centimeters','Position',[15,5,18,20]);
aw = 16;
ah = 2.8;
x0 = 1.5; dy = 1;
x = wirdum.y_coor/1000;
cmap = colormap('cool');

% Dip 
ax(1) = axes('Units','centimeters','Position',[x0, 4*ah + 5*dy, aw, ah]);
plot(x, input_table.dip); hold on
hwe = scatter(wirdum_eq.fault_rdy/1000, ones(size(wirdum_eq.fault_rdy)) *80, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Dip (deg)')

ax(2) = axes('Units','centimeters','Position',[x0, 3*ah + 4*dy, aw, ah]);
plot(x, input_table.dip_azi+90); hold on
plot([min(x),max(x)],[140,140]);
hwe = scatter(wirdum_eq.fault_rdy/1000, ones(size(wirdum_eq.fault_rdy)) *160, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Strike (deg)');     % this is dip azimuth, double check input of Vincent

ax(3) = axes('Units','centimeters','Position',[x0, 2*ah + 3*dy, aw, ah]);
plot(x, input_table.throw); hold on
hwe = scatter(wirdum_eq.fault_rdy/1000, ones(size(wirdum_eq.fault_rdy)) *70, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Throw (m)');   

ax(4) = axes('Units','centimeters','Position',[x0, 1*ah + 2*dy, aw, ah]);
plot(x, analysis_result.summary.reactivation_dp); hold on
hwe = scatter(wirdum_eq.fault_rdy/1000, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('\DeltaP_{reac} (MPa)');     

ax(5) = axes('Units','centimeters','Position',[x0,  1*dy, aw, ah]);
hold on
plot(x, analysis_result.summary.nucleation_dp);
hwe = scatter(wirdum_eq.fault_rdy/1000, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('\DeltaP_{nuc} (MPa)');     

xlabel(ax, 'Y RD (km)');
set(ax,'Box','on');
set([ax(4), ax(5)],'YLim',[-35,-8]);

%fig_path = '\\tsn.tno.nl\data\sv\sv-053185\Kluis\Research\DeepNL\PhysMax\Widrum_nucleation_modeling\';
%print(gcf,[fig_path, f_name, '.png'],'-dpng','-r300');

%% plot dip, strike, throw, reactivation and nucleation pressure along fault with events
% plotted against X RD
h2 = figure(2); clf(h2);
set(gcf, 'Units', 'centimeters','Position',[15,5,18,20]);
aw = 16;
ah = 2.8;
x0 = 1.5; dy = 1;
x = wirdum.x_coor/1000;
cmap = colormap('cool');

% Dip 
ax(1) = axes('Units','centimeters','Position',[x0, 4*ah + 5*dy, aw, ah]);
plot(x, input_table.dip); hold on
hwe = scatter(wirdum_eq.fault_rdx/1000, ones(size(wirdum_eq.rdx)) *80, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Dip (deg)')

ax(2) = axes('Units','centimeters','Position',[x0, 3*ah + 4*dy, aw, ah]);
plot(x, input_table.dip_azi+90); hold on
plot([min(x),max(x)],[140,140]);
hwe = scatter(wirdum_eq.fault_rdx/1000, ones(size(wirdum_eq.rdx)) *160, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Strike (deg)');     % this is dip azimuth, double check input of Vincent

ax(3) = axes('Units','centimeters','Position',[x0, 2*ah + 3*dy, aw, ah]);
plot(x, input_table.throw); hold on
hwe = scatter(wirdum_eq.fault_rdx/1000, ones(size(wirdum_eq.rdy)) *70, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('Throw (m)');   

ax(4) = axes('Units','centimeters','Position',[x0, 1*ah + 2*dy, aw, ah]);
plot(x, analysis_result.summary.reactivation_dp); hold on
hwe = scatter(wirdum_eq.fault_rdx/1000, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('\DeltaP_{reac} (MPa)');     

ax(5) = axes('Units','centimeters','Position',[x0,  1*dy, aw, ah]);
hold on
plot(x, analysis_result.summary.nucleation_dp);
hwe = scatter(wirdum_eq.fault_rdx/1000, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe,'MarkerEdgeColor','k');
ylabel('\DeltaP_{nuc} (MPa)');     

xlabel(ax, 'X RD (km)');
set(ax,'Box','on');
set([ax(4), ax(5)],'YLim',[-35,-22]);

fig_path = '\\tsn.tno.nl\data\sv\sv-053185\Kluis\Research\DeepNL\PhysMax\Widrum_nucleation_modeling\';
print(gcf,[fig_path, f_name, '_XRD_zoom.png'],'-dpng','-r300');


%% 

% mapview fault and events
h7 = figure(7); clf(h7); 
plot(wirdum.x_coor/1000, wirdum.y_coor/1000);
hold on
colormap cool
hwe1 = scatter(wirdum_eq.rdx/1000, wirdum_eq.rdy/1000, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
hwe2 = scatter(wirdum_eq.fault_rdx/1000, wirdum_eq.fault_rdy/1000, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled'); 
set(hwe2,'MarkerEdgeColor','k');
axis equal
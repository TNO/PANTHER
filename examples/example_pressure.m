% Example to explore how pressure is set with the Pressure object
% Pressure object sets initial pressures and pressures during load steps 
% in the FW, HW, and relative to those, in the fault
% In this example a figure is generated with FW HW and fault pressure
% Loes Buijze 05 - 01 - 2024

analysis = PantherInput();                          % intialize model run
analysis.input_parameters.throw.value = 50;
analysis.input_parameters.P_grad_res.value = 0.2;   % assign gas pressure gradient in the reservoir
analysis.input_parameters.P_over.value = 2;         % overpressure of 2 MPa at the top of the reservoir
analysis.input_parameters.width_FW.value = Inf;     % width of the footwall compartment
analysis.input_parameters.width_HW.value = Inf;     % width of the hanging wall compartment
analysis.P0_fault_mode = 'max';                     % initial fault pressure setting, relative to FW and HW pressures 
analysis.P_fault_mode = 'mean';                     % fault pressure during load steps, based on FW and HW pressures
analysis.load_table = analysis.load_table(1:3,:);   % reduce load table to 3 steps
analysis.load_table.time_steps(2) = 1;              % 1 year
analysis.load_table.P_steps(2) = -1;                % -1 MPa depletion
analysis.load_table.time_steps(3) = 30;             % 30 years
analysis.load_table.P_steps(3) = -30;               % -30 MPa depletion
analysis.diffusion_P = 1;                           % diffusion on (=1) or off (=0)

result = panther(analysis);
plot_pressures(result, analysis, height(analysis.load_table));

function [h2] = plot_pressures(result, analysis, t_step)

    h2 = figure(2); clf(h2); 
    set(gcf,'Units','centimeters','Position',[15,10,15,10]);
    clear ax
    aw = 4;
    ah = 8.5;
    dax = 0.5;
    ax0 = 1.5;
    ay0 = 1;
    % t_step = 2;
    run_description = ['\DeltaP: ', num2str(result.load_table.P_steps(t_step),'%.0f'),' MPa, ',...
        'w_{FW}: ', num2str(result.ensemble.width_FW(1),'%.0f' ), ' m, ', ...
        'w_{HW}: ', num2str(result.ensemble.width_HW(1),'%.0f' ), ' m, ', ...
        'p0 mode: ', analysis.P0_fault_mode,', ',...
        'p mode: ', analysis.P_fault_mode];

    figure_name = ['Fig_LargethrowDiff1_p0mode_',analysis.P0_fault_mode,'_pmode_',analysis.P_fault_mode,...
        '_dp_',num2str(result.load_table.P_steps(t_step),'%.0f'),...
        '_wFW_',num2str(result.ensemble.width_FW(1),'%.0f' ),...
        '_wHW_',num2str(result.ensemble.width_HW(1),'%.0f' )];

    P_HW = result.pressure{1}.get_P_HW();
    P_FW = result.pressure{1}.get_P_FW();
    
    dP_HW = result.pressure{1}.get_dP_HW();
    dP_FW = result.pressure{1}.get_dP_FW();
    
    % Initial pressures
    ax_num = 1; 
    
    ax(ax_num) = axes('Units','centimeters','Position',[ax0 + (ax_num-1)*(aw+dax),ay0, aw, ah]);
    plot(result.pressure{1}.P(:,1), result.y ,'Color','k','LineWidth', 1.5);
    hold on
    plot(P_HW(:,1), result.y , 'LineStyle','-.' ,'LineWidth', 1.5);
    plot(P_FW(:,1), result.y, 'LineStyle','--' ,'LineWidth', 1.5);
    grid on
    legend({'P0_{fault}','P0_{HW}','P0_{FW}'});
    xlabel('Initial pressure p0 (MPa)');
    

    % Absolute pressures
    ax_num = ax_num + 1;
    
    ax(ax_num) = axes('Units','centimeters','Position',[ax0 + (ax_num-1)*(aw+dax),ay0, aw, ah]);
    plot(result.pressure{1}.P(:,t_step), result.y ,'Color','k','LineWidth', 1.5);
    hold on
    plot(P_HW(:,t_step), result.y , 'LineStyle','-.' ,'LineWidth', 1.5);
    plot(P_FW(:,t_step), result.y, 'LineStyle','--' ,'LineWidth', 1.5);
    grid on
    legend({'P_{fault}','P_{HW}','P_{FW}'});
    xlabel('Pressure (MPa)');
    %title(['\DeltaP = ', num2str(result.load_table.P_steps(t_step),'%.0f'),' MPa']);
    title(run_description,'FontSize',8, 'HorizontalAlignment','center');

    % Pressure changes
    ax_num = ax_num + 1;
    ax(ax_num) = axes('Units','centimeters','Position',[ax0 + (ax_num-1)*(aw+dax),ay0, aw, ah]);
    plot(result.pressure{1}.dP(:,t_step), result.y ,'Color','k','LineWidth', 1.5);
    hold on
    plot(dP_HW(:,t_step), result.y , 'LineStyle','-.' ,'LineWidth', 1.5);
    plot(dP_FW(:,t_step), result.y, 'LineStyle','--' ,'LineWidth', 1.5);
    grid on
    legend({'\DeltaP_{fault}','\DeltaP_{HW}','\DeltaP_{FW}'});
    xlabel('Pressure change \Delta P (MPa)');
    
    linkaxes(ax,'y');
    set(ax(2:end),'YTickLabel','');
 
end

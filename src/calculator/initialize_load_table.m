function [load_table] = initialize_load_table()
    % initialize_load_table Initializes a load table with t, P, and T steps   
    % default 50 years, with 50 MPa depletion 
    % return: load_table
    % load_table.time_steps     [yr] time steps
    % load_table.P_steps        [MPa] pressure steps dP
    % load_table.P_factor_FW    [MPa] factor with which to multiply dP 
    % in the footwall compartment
    % load_table.P_factor_HW    [MPa] factor with which to multiply dP 
    % in the hanging wall compartment
    % load_table.P_factor_fault [MPa] factor with which to multiply dP 
    % in the fault
    % load_table.T_steps        [deg C] temperature steps dT
    % load_table.T_factor_FW    [MPa] factor with which to multiply dT 
    % in the footwall compartment
    % load_table.T_factor_HW    [MPa] factor with which to multiply dT 
    % in the hanging wall compartment

    load_table = table();
    load_table.time_steps = linspace(0, 35, 35)';
    n_steps = length(load_table.time_steps);
    load_table.P_steps = -1* linspace(0, 35, n_steps)';
    load_table.P_factor_HW = ones(size(load_table.P_steps));
    load_table.P_factor_FW = ones(size(load_table.P_steps));
    load_table.P_factor_fault = ones(size(load_table.P_steps));
    load_table.T_steps = -1* linspace(0, 35, n_steps)';
    load_table.T_multiply = ones(size(load_table.P_steps));
    load_table.T_factor_HW = ones(size(load_table.P_steps));
    load_table.T_factor_FW = ones(size(load_table.P_steps));

end
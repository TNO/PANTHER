function [run_results] = panther(analysis) 
    % main calculation. calculates fault stresses resulting from depletion
    % of one or two reservoir compartments, offset by a fault
    %
    % analysis: instance of input object generated by calling
    % PantherInput()

    % load default input parameters if no input provided
    if nargin == 0
        analysis = PantherInput();      
    end

    % generate the model ensemble, if it is not already generated
    if ~analysis.ensemble_generated 
        analysis.generate_ensemble;
    end

    % perform validation of input
    % TODO: validation functionality to be written

    % unpack some parameters
    y = analysis.y;
    load_table = analysis.load_table;
    load_case = analysis.load_case;
    p_fault = analysis.p_fault;
    p_res_mode = analysis.p_res_mode;
    diffusion_P = analysis.diffusion_P;
    diffusion_T = analysis.diffusion_T;
    ensemble = analysis.ensemble;
    dy = y(1) - y(2);
 
    % initialize output object
    run_results = PantherResult(y);

    % store table with input values in results
    run_results.ensemble = analysis.ensemble_to_table;

    % number of ensemble members
    n_members = length(analysis.ensemble);
    
    % initialize output arrays
    pressure = cell(n_members,1);
    temperature = cell(n_members,1);
    stress = cell(n_members,1);
    slip = cell(n_members,1);

    % set number of workers for parallel processing
    if or(n_members < 10, ~analysis.parallel)
        matlab_workers = 0;
    else
        matlab_workers = inf;
    end

    % calculate pressure and stress changes (parallel)
    % starting up the Matlab parallel pool (parpool) may take several 10s
    % of seconds.If it is already running (type gcp to check), parfor will initiate much faster
    parfor (i = 1 : n_members, matlab_workers)
        disp([num2str(i),'/', num2str(n_members)]);
        cell_length{i} = dy/sin(ensemble{i}.dip*pi/180);

        % initial stress
        initial_stress{i} = InitialStress(y, ensemble{i});
        
        % pressure and temperature changes
        pressure{i} = PantherPressure(ensemble{i}, y, load_table, load_case, diffusion_P, p_fault, p_res_mode);
        temperature{i} = Temperature(ensemble{i}, y, load_table, diffusion_T, 'min');
        
        % stress changes
        stress_change{i} = FaultStressChange(length(y), size(pressure{i}.dp_fault,2));        % initialize fault stresses for P
        stress_change{i} = stress_change{i}.calc_stress_changes(ensemble{i}, y, analysis.dx, pressure{i}, temperature{i}, load_case);
        
        % stress (initial + change)
        stress{i} = FaultStress(length(y), size(pressure{i}.dp_fault,2));
        stress{i} = stress{i}.compute_fault_stress(initial_stress{i}, stress_change{i}, pressure{i}.p);
        
        % fault slip, reactivation, nucleation
        slip{i} = FaultSlip(size(stress{i}.sne, 1), size(stress{i}.sne, 2));
        if analysis.aseismic_slip
            fault_strength{i} = stress{i}.sne.*ensemble{i}.f_s + ensemble{i}.cohesion;
            [slip{i}, stress{i}.tau] = slip{i}.calculate_fault_slip(y, stress{i}.sne, stress{i}.tau, ...
                                                         fault_strength{i}, ensemble{i}.get_mu_II);
        end
        
        slip{i} = slip{i}.detect_nucleation(cell_length{i}, stress{i}.sne, stress{i}.tau, ensemble{i}.f_s, ...
                                                    ensemble{i}.f_d, ensemble{i}.d_c, ensemble{i}.cohesion, ...
                                                    ensemble{i}.get_mu_II);
        % clear to save memory
        stress_change{i} = [];
        initial_stress{i}= [];
    end

    for i = 1 : n_members
        % save stresses, pressure, slip
        if analysis.save_stress
            run_results.pressure{i} = pressure{i};
            run_results.temperature{i} = temperature{i};
            run_results.stress{i} = stress{i};
            run_results.slip{i} = slip{i};
        end
    end

    % obtain summary report
    run_results = run_results.make_result_summary(analysis);


end
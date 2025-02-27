function [analysis] = panther(analysis) 
    % main calculation. calculates fault stresses resulting from depletion
    % of one or two reservoir compartments, offset by a fault
    % can handle stochastic calculation with all input parameters in the
    % ensemble member, but not for depth-dependent properties
    % all calculations have the same run settings (e.g. y, diffusion_P)
    % to run simulations for faults with different depth-dependent
    % properties or different run settings, use MultiFaultCalculator
    % INPUT
    % analysis: instance of input object generated by calling
    % PantherInput()

    % load default input parameters if no input provided
    if nargin == 0
        analysis = PantherInput();      
    end

    % (re)generate ensemble member objects
    analysis.generate_ensemble();

    % validate input
    validate_input(analysis);

    % unpack some parameters
    y = analysis.y;
    load_table = analysis.load_table;
    load_case = analysis.load_case;
    diffusion_T = analysis.diffusion_T;
    ensemble_members = analysis.ensemble_members;
    nucleation_criterion = analysis.nucleation_criterion;
    nucleation_length = analysis.nucleation_length_fixed;
    suppress_status_output = analysis.suppress_status_output;

    % define output steps
    if ~ismember(analysis.save_stress,'none')
           if ismember(analysis.save_stress,'all')
               indices_for_saving = 1 : height(load_table);
           elseif ismember(analysis.save_stress,'first')
               indices_for_saving = 1;
           elseif ismember(analysis.save_stress,'last')
                indices_for_saving = size(height(load_table),2);
           elseif ismember(analysis.save_stress,'first_last')
               indices_for_saving = [1, size(height(load_table),2)];
           end
    else
        indices_for_saving = [];
    end

    % number of ensemble members
    n_members = length(analysis.ensemble_members);
    
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
    % starting up the Matlab parallel pool (parpool) for the first time may take several 10s
    % of seconds. If it is already running (type gcp to check), parfor will initiate much faster
    parfor (i = 1 : n_members, matlab_workers)
        if ~suppress_status_output
            disp([num2str(i),'/', num2str(n_members)]);
        end
        L{i} = y./sin(ensemble_members{i}.dip*pi/180);

        % initial stress
        initial_stress{i} = InitialStress(y, ensemble_members{i});
        
        % pressure and temperature changes
        pressure{i} = Pressure(ensemble_members{i}, load_table, analysis);
        temperature{i} = Temperature(ensemble_members{i}, y, load_table, diffusion_T, 'min');
        
        % stress changes
        stress_change{i} = FaultStressChange(length(y), size(pressure{i}.dP,2));        % initialize fault stresses for P
        stress_change{i} = stress_change{i}.calc_stress_changes(ensemble_members{i}, y, analysis.dx, pressure{i}, temperature{i}, load_case);
        
        % stress (initial + change)
        stress{i} = FaultStress(length(y), size(pressure{i}.dP,2));
        stress{i} = stress{i}.compute_fault_stress(initial_stress{i}, stress_change{i}, pressure{i}.P);
        
        % fault slip, reactivation, nucleation
        slip{i} = FaultSlip(size(stress{i}.sne, 1), size(stress{i}.sne, 2));
        if analysis.aseismic_slip
            fault_strength{i} = stress{i}.sne.*ensemble_members{i}.f_s + ensemble_members{i}.cohesion;
            [slip{i}, stress{i}.tau] = slip{i}.calculate_fault_slip(L{i}, stress{i}.sne, stress{i}.tau, ...
                                                         fault_strength{i}, ensemble_members{i}.get_mu_II);
        end
        
        slip{i} = slip{i}.detect_nucleation(L{i}, stress{i}.sne, stress{i}.tau, ensemble_members{i}.f_s, ...
                                                    ensemble_members{i}.f_d, ensemble_members{i}.d_c, ensemble_members{i}.cohesion, ...
                                                    ensemble_members{i}.get_mu_II, nucleation_criterion, nucleation_length);
        % clear to save memory
        stress_change{i} = [];
        initial_stress{i}= [];

        % get the fault stressses at onset of reactivation and nucleation 
        stress{i} = stress{i}.get_reactivation_stress(slip{i}.reactivation_load_step);
        stress{i} = stress{i}.get_nucleation_stress(slip{i}.nucleation_load_step);
        
        % reduce output
        pressure{i} = pressure{i}.reduce_steps(indices_for_saving);
        stress{i} = stress{i}.reduce_steps(indices_for_saving);
        temperature{i} = temperature{i}.reduce_steps(indices_for_saving);
        slip{i} = slip{i}.reduce_steps(indices_for_saving);

    end
    
    reduced_load_table = analysis.load_table(indices_for_saving,:);   

    for i = 1 : n_members
        % save stresses, pressure, slip
        if ~ismember(analysis.save_stress,'none')
            analysis.pressure{i} = pressure{i};
            analysis.temperature{i} = temperature{i};
            analysis.stress{i} = stress{i};
            analysis.slip{i} = slip{i};
            analysis.load_table = reduced_load_table;
        end
    end

    % obtain summary report
    analysis = analysis.make_result_summary();
    
end
function validate_input(analysis)
    % function to validate input which has not yet been validated in the
    % class properties, and perform sanity checks on the chosen input
    % parameters
    
    % check that the correct input is entered for the save_stress option
    if ischar(analysis.save_stress{1})
        if ~ismember(analysis.save_stress{1},{'all','none','first_last','first','last'});
             error(['For save stress selection enter step numbers or ',...
                'all, none, first_last, first, last' ]);
        end
    end
    
    % check sensibility of the chosen Young's modulus
    if log10(analysis.input_parameters.young.value) < 3 || log10(analysis.input_parameters.young.value) > 5
        warning('NB default input unit Youngs modulus is MPa, please double check your input');
    end

    % check sensibility of the chosen vertical stress gradient
    if analysis.input_parameters.sv_grad.value < 18 || analysis.input_parameters.sv_grad.value > 28
        warning(['Vertical stress gradient in in MPa/km, and lies typically ',...
            'in the range of 19 - 24 MPa/km. Given input was ', ...
            num2str(analysis.input_parameters.sv_grad.value),' please check']);
    end

    % check sensibility of the chosen pressure gradient
    if analysis.input_parameters.P_grad.value < 9 || analysis.input_parameters.P_grad.value > 15
        warning(['Pore stress gradient in in MPa/km, and lies typically ',...
            'in the range of 10 - 12 MPa/km. Given input was ', ...
            num2str(analysis.input_parameters.P_grad.value),' please check']);
    end

    % check sensibility of the chosen temperature gradient
    if analysis.input_parameters.T_grad.value < 26 || analysis.input_parameters.T_grad.value > 40
        warning(['Pore stress gradient in in deg/km, and lies typically ',...
            'in the range of 26 - 35 deg/km far away from plate boundaries.',...
            'Given input was ', ...
            num2str(analysis.input_parameters.T_grad.value),' please check']);
    end

    % check sensibility of the chosen depth
    if analysis.input_parameters.depth_mid.value < -5000 
        warning(['Depth of subsurface activities is typically < 5000 m deep ',...
            'Given input was ', ...
            num2str(analysis.input_parameters.depth_mid.value),' please check']);
    end

    % check for error in throw
    if analysis.input_parameters.throw.value > abs(analysis.input_parameters.depth_mid.value)
        error(['Throw of ', num2str(analysis.input_parameters.throw.value),...
            ' exceeds reservoir depth']);
    end

    % check for error in throw
    if analysis.input_parameters.thick.value > abs(analysis.input_parameters.depth_mid.value)
        error(['Thickness of ', num2str(analysis.input_parameters.thick.value),...
            ' exceeds reservoir depth']);
    end

    depth_dependent_properties = analysis.input_parameters.get_depth_dependent_properties();
    allowed_depth_dependent_properties = {'dip','dip_azi','young','poisson',...
        'biot','therm_exp','sH_dir','shsv','sHsh','f_s','f_d','d_c','cohesion'};
        if ~isempty(depth_dependent_properties)
            faulty_input = ~ismember(depth_dependent_properties, allowed_depth_dependent_properties);
            if any(faulty_input)
                i_first_faulty = find(faulty_input, 1, 'first');
                error(['Depth dependent property set where not allowed ',...
                    'e.g. ', depth_dependent_properties{i_first_faulty}, ...
                    '. Depth dependent properties are dip,dip_azi, young, poisson,',...
                    'biot,therm_exp,sH_dir,shsv,sHsh,f_s, f_d,d_c,cohesion']);
            end
            % set stochastic mode for depth dependent properties off
            for i = 1 : length(depth_dependent_properties)
                input_value_length = length(analysis.input_parameters.(depth_dependent_properties{i}).value_with_depth);
                if (length(analysis.y) ~= length(input_value_length)) & (~input_value_length == 1 )
                    error(['Length of depth dependent property ', depth_dependent_properties{i},...
                        ' does not seem to match length of y, please check input']);
                elseif input_value_length == 1 
                    warning(['Depth dependency set for ', depth_dependent_properties{i},...
                        ' but length of value with depth is 1']);
                end
                analysis.input_parameters.(depth_dependent_properties{i}).stochastic = 0;
            end
        end

end
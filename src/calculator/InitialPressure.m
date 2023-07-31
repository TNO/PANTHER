classdef InitialPressure
    % InitialPressure sets the initial pore pressure in the reservoir
    % compartments, seal and base, and in the fault. Within the reservoir
    % the fluid can have a different density. 

    properties
        p0 double       % [MPa] initial pressure fault
        p0_HW double    % [MPa] initial pressure hanging wall compartment
        p0_FW double    % [MPa] initial pressure footwall compartment
    end

    methods
        function self = InitialPressure(member, y, p_fault_mode, p_res_mode)
            [next_to_FW, next_to_HW, ~] = is_adjacent_to_reservoir(y, member.thick, member.throw); 
            yy = y + member.depth_mid;    
            p0 = zeros(size(yy));                       % intialize fault pressure
            p0 = -(yy/1000).*member.p_grad + member.p_offset;  % [MPa] hydrostatic pressure
            p0_FW = p0; p0_HW = p0;                     % initialize FW and HW pressure with hydrostatic
            top_HW_i = find(next_to_HW,1, 'first');     % index where HW compartment starts (top)
            %top_HW_i = member.top_HW_i(y)
            top_FW_i = find(next_to_FW,1, 'first');     % index where FW compartment starts (top)
            top_res_i = min(top_HW_i, top_FW_i);        % top most depth of reservoir interval
            base_FW_i = find(next_to_FW,1, 'last'); 
            base_HW_i = find(next_to_HW,1, 'last');
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base
            p0_HW(top_HW_i:end) = p0_HW(top_HW_i:end) + member.p_over;  
            p0_FW(top_FW_i:end) = p0_FW(top_FW_i:end) + member.p_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(p_res_mode, 'same')
                base_res_i = max(base_FW_i, base_HW_i );
                p0_FW(top_FW_i:base_res_i) = p0_FW(base_res_i) - (yy(top_FW_i:base_res_i) - yy(base_res_i))*member.p_grad_res/1000;
                p0_HW(top_HW_i:base_res_i) = p0_HW(base_res_i) - (yy(top_HW_i:base_res_i) - yy(base_res_i))*member.p_grad_res/1000;
            else
                p0_FW(next_to_FW) = p0_FW(base_FW_i) - (yy(next_to_FW) - yy(base_FW_i))*member.p_grad_res/1000;
                p0_HW(next_to_HW) = p0_HW(base_HW_i) - (yy(next_to_HW) - yy(base_HW_i))*member.p_grad_res/1000; 
            end
            % set fault pressure
            if strcmp(p_fault_mode, 'max')
                p0 = max(p0_HW, p0_FW);
            elseif strcmp(p_fault_mode, 'min')
                p0 = min(p0_HW, p0_FW);
                p0 = max(p0_HW, p0_FW); % TODO fix bug here, for initial pressure it should be the max p 
            elseif strcmp(p_fault_mode, 'mean')
                p0 = mean([p0_HW, p0_FW],2);
            elseif strcmp(p_fault_mode, 'FW')
                p0 = p0_FW;
            elseif strcmp(p_fault_mode, 'HW')
                p0 = p0_HW;
            end
            self.p0 = p0;
            self.p0_HW = p0_HW;
            self.p0_FW = p0_FW;
        end
    end

end
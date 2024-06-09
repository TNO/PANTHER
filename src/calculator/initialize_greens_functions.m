function [GF] = initialize_greens_functions(params, y, dx, variable_PT, variable_dip)
    % InitializeGreens Initializes Green's functions for run
    % INPUT
    % params                input parameters ensemble member
    % y                     depth w.r.t. y_mid
    % dx                    distance from fault in x
    % varies_with_depth     1 if P and / or T vary with y
    % depth
    % OUTPUT
    % GF    cell array of one Green's function (in case of uniform) or of
    %       size (len(y), 1) in case of varying P,T,dip,or elastic properties
    %       with depth

  
    % numerical_correction Corrects x and y in case these coincide
    % with any of the boundaries of the reservoir shapes. this
    % results in singularities and division by 0
    correction_value = 1e-9;
    if dx == 0 || dx == params.width_FW || dx == -params.width_HW
        dx = dx + correction_value;
    end
    reservoir_boundaries = [params.top_FW_y, params.top_HW_y, params.base_HW_y, params.base_FW_y];
    if any(ismember(y, reservoir_boundaries))
        i_boundary = find(ismember(y, reservoir_boundaries));
        y(i_boundary) = y(i_boundary) + correction_value;
    end
    % xcoordinates
    xeval = y./(tan(params.dip*pi/180)) + dx;
        if and(~variable_PT, ~variable_dip)
            GF{1} = GreensFunctions(y);
            if params.width_FW > 0  % and add a criterium to check if the FW dP and dT are not 0
                GF{1} = GF{1}.green_FW(xeval, y, params.dip, params.thick, params.throw, params.width_FW, 0, 0 );
            end
            if params.width_HW > 0 % and add a criterium to check if the FW dP and dT are not 0
                GF{1} = GF{1}.green_HW(xeval, y, params.dip, params.thick, params.throw, params.width_HW, 0, 0 );
            end
        elseif and(variable_PT, ~variable_dip)
            % calculate the same GF once, but shift it along depth axis
            slice_thick = (y(1) - y(2));                        % depth slice thickness
            y2 = [y; y(1:end-1)+(y(end)-y(1))-slice_thick];     % pad with zeros
            xeval2 = y2/(tan(params.dip*pi/180)) + dx;    
            i_mid = ceil(length(y2)/2);
            slice_y = y2(i_mid);                                % depth slice mid y on fault
            slice_x = slice_y/(tan(params.dip*pi/180));         % depth slice mid x on fault
            slice_throw = 0;                                    % depth slice throw, set to 0. 
            greens_f = GreensFunctions(y2);                % initialize Green's functions
            if params.width_FW > 0  % and add a criterium to check if the FW dP and dT are not 0
                greens_f = greens_f.green_FW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
            end
            if params.width_HW > 0 % and add a criterium to check if the FW dP and dT are not 0
                greens_f = greens_f.green_HW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_HW, slice_x, slice_y);
            end
            GF = cell(length(y), 1);
            for j = 1 : length(y)
                GF{j} = greens_f;
                GF{j}.Gnorm_FW = GF{j}.Gnorm_FW(i_mid-j+1:2*i_mid-j);
                GF{j}.Gnorm_HW = GF{j}.Gnorm_HW(i_mid-j+1:2*i_mid-j);
                GF{j}.Gshear_FW = GF{j}.Gshear_FW(i_mid-j+1:2*i_mid-j);
                GF{j}.Gshear_HW = GF{j}.Gshear_HW(i_mid-j+1:2*i_mid-j);
                % here, remove the xx components if needed
            end
        else
            % calculate separate GF at each depth interval (slowest)
            % this is needed if the geometry (e.g. dip, w_HW etc) changes
            slice_thick = y(1) - y(2);                      % depth slice thickness
            slice_throw = 0; 
            GF = cell(length(y), 1);
            % depth slice throw, set to 0. 
            for j = 1 : length(y)
                if variable_dip
                    dip = params.dip(j);
                else
                    dip = params.dip;
                end
                slice_y = y(j); % - 0.5*slice_thick;        % depth slice mid y on fault
                slice_x = slice_y/(tan(dip*pi/180)); % depth slice mid x on fault
                GF{j} = GreensFunctions(y);              % initialize Green's functions
                GF{j} = GF{j}.green_FW(xeval, y, dip, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
                GF{j} = GF{j}.green_HW(xeval, y, dip,  slice_thick, slice_throw, params.width_HW, slice_x, slice_y );
            end
        end            

    end

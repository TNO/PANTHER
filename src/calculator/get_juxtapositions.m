% creates a table with FW, HW and combined FWHW stratigraphies
function [jux] = get_juxtapositions(varin)

    if ~isfield(varin,'TBthick')
        varin.TBthick = 50; % default 50 m Ten Boer thickness
    end
    ZEthick = 50;
    
    thickness = varin.thick;
    offset = varin.throw;
    
    a = (thickness - offset)/2;  % +ve top left compartment, -ve base right
    b = (thickness + offset)/2;  % 
    topFW = b; baseFW = -a;
    topHW = a; baseHW = -b;
    
    y = varin.y;        % depth value w.r.t. center depth, deep to shallow
    
    jux = table();
    jux.FW = cell(size(y));     % footwall stratigraphies
    jux.FW(y < baseFW) = {'DC'};                        % Carboniferous
    jux.FW(y >= baseFW & y < topFW) = {'SS'};           % Slochteren reservoir
    jux.FW(y >= topFW & y < (topFW + varin.TBthick)) = {'TB'}; % Ten Boer
    jux.FW(y >= (topFW + varin.TBthick) & y < (topFW + varin.TBthick + ZEthick)) = {'BZ'};
    jux.FW(y >= (topFW + varin.TBthick + ZEthick)) = {'ZE'}; % Zechstein
    jux.HW = cell(size(y));     % hangin wall stratigraphies
    jux.HW(y < baseHW) = {'DC'};
    jux.HW(y >= baseHW & y < topHW) = {'SS'};
    jux.HW(y >= topHW & y < (topHW + varin.TBthick)) = {'TB'};
    jux.HW(y >= (topHW + varin.TBthick) & y < (topHW + varin.TBthick + ZEthick)) = {'BZ'};
    jux.HW(y >= (topHW + varin.TBthick + ZEthick)) = {'ZE'};
    
    jux.FWHW = strcat(jux.FW, jux.HW);

end
classdef LoadFigure

    properties
        figure_width (1,1) double {mustBeNonnegative} = 16
        figure_height (1,1) double {mustBeNonnegative} = 8.5
        axes_width (1,1) double {mustBeNonnegative} = 3
        axes_height (1,1) double {mustBeNonnegative} = 7
        axes_x0 (1,1) double {mustBeNonnegative} = 1.5
        axes_y0 (1,1) double {mustBeNonnegative} = 1
        axes_xspacing (1,1) double {mustBeNonnegative} = 0.5
        axes_yspacing (1,1) double {mustBeNonnegative} = 1
        axes_setting {mustBeMember(axes_setting, {'explicit','auto'})} = 'explicit';
        n_subplots {mustBeInteger, mustBeNonnegative} = 4
        n_rows {mustBeInteger, mustBeNonnegative} = 1
        n_columns {mustBeInteger, mustBeNonnegative} = 4
        units {mustBeMember(units, {'inches','pixels','centimeters', 'normalized','points','characters'})} = 'centimeters'
        annotate logical = true
        annotation_type {mustBeMember(annotation_type, {'alphabetic','numeric','roman'})} = 'alphabetic'
        annotation_unit {mustBeMember(annotation_unit, {'inches','pixels','centimeters', 'normalized','points','characters'})} = 'normalized'
        annotation_x = 0.85
        annotation_y = 0.05
        annotation_alignment {mustBeMember(annotation_alignment,{'center','left','right'})} = 'center'
        axes_font_size double {mustBeNonnegative} = 9
    end

    methods

        function self = LoadFigure()
        end

        function [h] = load(self)
            % result: result object
            h_open_figs =  findobj('type','figure');
            number_of_open_figures = length(h_open_figs);
            % open figures
            h.fig = figure(number_of_open_figures + 1); clf(h.fig);
            set(gcf, 'Units', self.units, 'Position',[15,10, self.figure_width, self.figure_height]);
            movegui(h.fig, 'center');
            h.ax = self.set_subplot_axes();
            if self.annotate
                h.annotations = self.annotate_subplots(h.ax);
            end
            set(h.ax,'FontSize', self.axes_font_size);
        end

        function [ax_handle] = set_subplot_axes(self)
            % Create number of subplot axes, specified by n_rows and n_columns. 
            x0 = self.axes_x0;          % axes bottom left corner x
            y0 = self.axes_y0;          % axes bottom left corner y
            ax_w = self.axes_width;        % axes width
            ax_h = self.axes_height;       % axes height
            dx = self.axes_xspacing;    % axes spacing x
            dy = self.axes_yspacing;    % axes spacing y
        
            fig_width = self.figure_width;
            fig_height = self.figure_height;
        
            % automatically rescale axes if exceeding figure size
            if (self.n_rows*ax_h + (self.n_rows-1)*dy + y0) > fig_height
                y0 = 0.2 * (fig_height/self.n_rows);
                ax_h = 0.75 * (fig_height/self.n_rows);
                dy = 0.2 * (fig_height/self.n_rows);
                disp('NB Axes plotting outside figure bounds, adjust axes heigth');
            end
            if (self.n_columns*ax_w + (self.n_columns-1)*dx + x0) > fig_width
                x0 = 0.15 * (fig_width/self.n_columns);
                ax_w = 0.8 * (fig_width/self.n_columns);
                dx = 0.15 * (fig_width/self.n_columns);
                disp('NB Axes plotting outside figure bounds, adjust axes width');
            end
                  
            for i = 1 : self.n_rows*self.n_columns
                row = ceil(i/self.n_columns);   % row_number
                col = i - self.n_columns*(row-1);
                x0_ax = x0 + (col-1) * (dx + ax_w);
                y0_ax = y0 + (self.n_rows - row)*(dy + ax_h);
                ax_handle(i) = axes('Units','centimeters','Position',[x0_ax, y0_ax, ax_w, ax_h],'Box','on',...
                    'FontSize',self.axes_font_size);
            end
        end

        function [axes_handles] = set_axes_auto(self)
            % calculate axes width and height automatically

        end

        function annotation_handles = annotate_subplots(self, axes_handles)
           annotations = self.load_annotation_array();
            for i = 1 : length(axes_handles)
                axes(axes_handles(i));
                annotation_handles(i) = text(self.annotation_x, self.annotation_y, annotations{i}, 'Units',self.annotation_unit);
            end
            set(annotation_handles,'HorizontalAlignment', self.annotation_alignment);
        end

        function annotation_array = load_annotation_array(self)
             if strcmp(self.annotation_type,'alphabetic')
                annotation_array = {'a', 'b','c','d','e','f','g','h','i','j','k','l','m',...
                'n','o','p','q','r','s','t','u','v','w','x','y','z'};
             elseif strcmp(self.annotation_type,'numeric')
                 numeric_array = 1:1:100;
                 annotation_array = cellstr(num2str(numeric_array(:)))';
             elseif strcmp(self.annotation_type,'roman')
                numeric_array = 1:1:100;
                annotation_array = self.num2roman(numeric_array);
             end
        end

        function romanNumerals = num2roman(~, numeric_array)
            % Define the Roman numeral symbols and their corresponding values
            romanSymbols = {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
            romanValues = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
            
            % Initialize the output cell array
            romanNumerals = cell(size(numeric_array));
            
            % Loop through each number in the array
            for i = 1:length(numeric_array)
                num = numeric_array(i);
                romanStr = '';
                
                % Convert the number to Roman numeral
                for j = 1:length(romanValues)
                    while num >= romanValues(j)
                        romanStr = [romanStr, romanSymbols{j}]; % Append the Roman symbol
                        num = num - romanValues(j);           % Subtract the corresponding value
                    end
                end
                
                romanNumerals{i} = romanStr; % Store the result in the cell array
            end
        end


    end

    
end
classdef Plot1DResult < LoadFigure
    % plots results on the 1D fault
    % does not yet work for stochastic, picks the first model run
    
    properties
        plot_results cell = {'P','sne','tau','scu','slip'};
        plot_initial logical = true
        plot_step = 'last'
        plot_colors = {'k'}
        plot_colors_initial = [0.5,0.5,0.5]
        plot_handles
        % side_indicator logical = true
        plot_linewidth double {mustBePositive} = 1.5
        plot_linestyle_initial = '--'
        ax_scale {mustBeMember(ax_scale,{'auto','explicit'})} = 'auto';
        xlim 
        ylim 
        grid_on logical = true
    end

    properties (Constant)
        
    end

    methods

        function self = Plot1DResult()
            % Initialize
            disp('Initialize plot for on-fault results');
        end

        function self = plot_PANTHER_result(self, analysis)
            % here
            self.n_rows = 1;
            self.n_columns = length(self.plot_results);
            self.axes_width = 2.5;
            self.axes_xspacing = 0.3;
            [self.plot_handles] = self.load();  % load default figure
            n_steps = analysis.nTimes;
            if strcmp(self.plot_step,'last')
                i_step = n_steps;
            elseif isnumeric(self.plot_step) && self.plot_step <= n_steps
                i_step = self.plot_step;
            else
                i_step = n_steps;
            end
            y_abs  = analysis.getDepth();
            for i = 1 : length(self.plot_results)
                x{i} = self.retrieve_result_plot_data(analysis, self.plot_results{i}, i_step);
                if self.plot_initial
                    x_ini{i} = self.retrieve_result_plot_data(analysis,  self.plot_results{i}, 1);
                end
            end
            x_labels = self.retrieve_labels();
            % loop over axes
            for i = 1 : length(self.plot_results)
                axes(self.plot_handles.ax(i));
                hold on
                % plot results
                if self.plot_initial
                    self.plot_handles.axis_plots_ini(i) = plot(x_ini{i}, y_abs);                
                end
                self.plot_handles.axis_plots(i) = plot(x{i}, y_abs); 
                % set colors. black if default if no color specified
                if length(self.plot_colors) >= i
                    set(self.plot_handles.axis_plots,'Color',self.plot_colors{i});
                else 
                    set(self.plot_handles.axis_plots,'Color','k');
                end
                % set axes labels
                xlabel(x_labels{i});
                if i == 1
                    ylabel('Depth (m)');
                end
                % set grid
                if self.grid_on
                    grid on
                end
                % set axes limits
                if strcmp(self.ax_scale,'auto')
                    set(gca, 'YLim', [min(y_abs), max(y_abs)]);
                else
                    if ~isempty(self.ylim) && strcmp(self.ax_scale, 'explicit')
                        set(gca,'YLim', self.ylim);
                    end
                    if ~isempty(self.xlim)
                        set(gca,'XLim', self.xlim);
                    end
                end
            end
            if length(self.plot_handles.ax) > 1
                set(self.plot_handles.ax(2:end), 'YTick','');
            end
            set(self.plot_handles.axis_plots,'LineWidth',self.plot_linewidth);
            if self.plot_initial
                set(self.plot_handles.axis_plots_ini,'LineWidth',self.plot_linewidth);
                set(self.plot_handles.axis_plots_ini,'Color',self.plot_colors_initial);
                set(self.plot_handles.axis_plots_ini,'LineStyle',self.plot_linestyle_initial);
            end
        end

        function [array_to_plot] = retrieve_result_plot_data(~, analysis, parameter, i_step)
            
            if contains(parameter, 'scu' )
                scu = analysis.getSCU();
                array_to_plot = scu(:, i_step);
            else
                array_to_plot = analysis.faultResults.(parameter)(:,i_step);
            end
        end

        function x_labels = retrieve_labels(self)
            x_labels = cell(size(self.plot_results));
            x_labels(:) = {''};
            for i = 1 : length(self.plot_results)
                if strcmp(self.plot_results{i}, 'P' )
                    x_labels{i} = 'Pressure (MPa)';
                elseif contains(self.plot_results{i}, 'sne' )
                    x_labels{i} = '\sigma_n'' (MPa)';
                elseif contains(self.plot_results{i}, 'tau' )
                    x_labels{i} = '\tau (MPa)';
                elseif strcmp(self.plot_results{i}, 'scu' )
                    x_labels{i} = 'SCU (-)';
                elseif strcmp(self.plot_results{i}, 'slip' )
                    x_labels{i} = 'Slip (m)';
                end
            end
        end
        
          function plot_side_indicators(~, ~)
        end

          function self = plot_custom_plot_in_axes(self)
        end

          function validate_plot_variables(~)
           % {mustBeMember(plot_results, 'P','T','sne','tau','sne_reac','tau_reac','sne_nuc','tau_nuc')}
        end
    end

end
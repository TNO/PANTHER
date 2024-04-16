function plot_PANTHER_results(input_table, wirdum, wirdum_eq, analysis_result, output_dir, direction)
    % Plot dip, strike, throw, reactivation, and nucleation pressure along fault with events
    h2 = figure(2); clf(h2);
    set(gcf, 'Units', 'centimeters','Position',[15,5,18,20]);
    aw = 16;
    ah = 2.8;
    x0 = 1.5; dy = 1;

    if direction == 'x'
        x_data = wirdum.x_coor/1000;
        xlabel_label = 'X RD (km)';
        scatter_x = wirdum_eq.fault_rdx/1000;
        scatter_y = wirdum_eq.rdx;
    elseif direction == 'y'
        x_data = wirdum.y_coor/1000;
        xlabel_label = 'Y RD (km)';
        scatter_x = wirdum_eq.fault_rdy/1000;
        scatter_y = wirdum_eq.rdy;
    else
        error('Invalid direction. Use either ''x'' or ''y''.');
    end

    cmap = colormap('cool');

    for i = 1:5
        ax(i) = axes('Units','centimeters','Position',[x0, (5-i)*ah + (6-i)*dy, aw, ah]);
        switch i
            case 1
                plot(x_data, input_table.dip); hold on
                ylabel('Dip (deg)');
                scatter(scatter_x, ones(size(scatter_x)) * 80, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled');
            case 2
                plot(x_data, input_table.dip_azi+90); hold on
                plot([min(x_data),max(x_data)],[140,140]);
                ylabel('Strike (deg)');
                scatter(scatter_x, ones(size(scatter_x)) * 160, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled');
            case 3
                plot(x_data, input_table.throw); hold on
                ylabel('Throw (m)');
                scatter(scatter_x, ones(size(scatter_x)) * 70, 10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled');
            case 4
                plot(x_data, analysis_result.summary.reactivation_dp); hold on
                ylabel('\DeltaP_{reac} (MPa)');
                scatter(scatter_x, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled');
            case 5
                plot(x_data, analysis_result.summary.nucleation_dp); hold on
                ylabel('\DeltaP_{nuc} (MPa)');
                scatter(scatter_x, wirdum_eq.pressure_change_at_eq_time,10.^(1+wirdum_eq.magnitude/2.2), wirdum_eq.years_from_first_event,'filled');
        end
        xlabel(xlabel_label);
    end

    set(ax,'Box','on');
    % set([ax(4), ax(5)],'YLim',[-35,-22]);

    % Save figure to the output directory
    filename = fullfile([output_dir direction '.png']);
    saveas(h2, filename);
end

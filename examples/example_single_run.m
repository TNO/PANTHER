% example file for a single model run

% initialize input 
run_instance = PantherInput;
% change some input properties
run_instance.input_parameters.width_HW.value = 500;
run_instance.input_parameters.throw.value = 50;
% turn off diffusion
run_instance.diffusion_P = 0;
% generate model ensemble to check input
run_instance.generate_ensemble();
% run panther with current input instance
run_instance = panther(run_instance);


hfig = Plot1DResult();
hfig.axes_font_size = 8;
hfig.ax_scale = 'explicit';
hfig.ylim = [run_instance.ensemble.depth_mid(1) - 300, run_instance.ensemble.depth_mid(1) + 300];
hfig.plot_PANTHER_result(run_instance.input_parameters, run_instance);



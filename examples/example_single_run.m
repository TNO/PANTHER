% example file for a single model run

% initialize input 
run_instance = PantherAnalysis();
% change some input properties
run_instance.setInputParameter('width_HW', 500);
run_instance.setInputParameter('throw', 50);
% turn off diffusion
run_instance.diffusion_P = 0;
% generate model ensemble to check input
run_instance.generate_ensemble();
% run panther with current input instance
run_instance.run();

hfig = Plot1DResult();
hfig.axes_font_size = 8;
hfig.ax_scale = 'explicit';
depth_mid = run_instance.getInputParameter('depth_mid');
hfig.ylim = [depth_mid - 300, depth_mid + 300];
hfig.plot_PANTHER_result(run_instance);



% example stochastic run

analysis = PantherInput;
analysis.input_parameters.dip.stochastic = 1;   % make dip a stochastic parameter
analysis.input_parameters.shsv.stochastic = 1;  % make shsv a stochastic parameter
analysis.input_parameters.shsv.a = 0.69;        % lower value of uniform distribution
analysis.input_parameters.shsv.b = 0.76;        % upper value of uniform distribution
analysis.stochastic = 1;                        % set the analysis to stochastic    
analysis.n_stochastic = 100;                     % number of stochastic runs
analysis.generate_ensemble();                   % generate the ensemble (run_instance.ensemble)

analysis_result = panther(analysis);


%%

analysis_input = analysis.ensemble_to_table();   

h1 = figure(1); clf(h1);
subplot(2,2,1)
scatter(analysis_input.dip, analysis_result.summary.nucleation_dp);
set(gca,'YDir','reverse','Box','on');
ylabel('Nucleation depletion (MPa)');
xlabel(['Dip (', char(176), ')'])
subplot(2,2,2)
scatter(analysis_input.shsv, analysis_result.summary.nucleation_dp);
set(gca,'YDir','reverse','Box','on');
ylabel('Nucleation depletion (MPa)');
xlabel('\sigma_h / \sigma_v (-)');
subplot(2,2,3)
scatter(analysis_input.dip, analysis_result.summary.mid_cff_rate);
set(gca,'Box','on');
ylabel('CFF rate mid (MPa/yr)');
xlabel(['Dip (', char(176), ')'])
subplot(2,2,4)
scatter(analysis_input.shsv, analysis_result.summary.mid_cff_rate);
set(gca,'Box','on');
ylabel('CFF rate mid (MPa/yr)');
xlabel('\sigma_h / \sigma_v (-)');
% example stochastic run

analysis = PantherAnalysis;
analysis.input_parameters.dip.stochastic = 1;   % make dip a stochastic parameter
analysis.input_parameters.shsv.stochastic = 1;  % make shsv a stochastic parameter
analysis.input_parameters.shsv.a = 0.69;        % lower value of uniform distribution
analysis.input_parameters.shsv.b = 0.76;        % upper value of uniform distribution
analysis.stochastic = 1;                        % set the analysis to stochastic    
analysis.n_stochastic = 20;                     % number of stochastic runs
analysis.generate_ensemble();
analysis = panther(analysis);


%%

h1 = figure(1); clf(h1);
subplot(2,2,1)
scatter(analysis.ensemble.dip, analysis.summary.nucleation_dP);
set(gca,'YDir','reverse','Box','on');
ylabel('Nucleation depletion (MPa)');
xlabel(['Dip (', char(176), ')'])
subplot(2,2,2)
scatter(analysis.ensemble.shsv, analysis.summary.nucleation_dP);
set(gca,'YDir','reverse','Box','on');
ylabel('Nucleation depletion (MPa)');
xlabel('\sigma_h / \sigma_v (-)');
subplot(2,2,3)
scatter(analysis.ensemble.dip, analysis.summary.cff_ymid);
set(gca,'Box','on');
ylabel('CFF rate mid (MPa/yr)');
xlabel(['Dip (', char(176), ')'])
subplot(2,2,4)
scatter(analysis.ensemble.shsv, analysis.summary.cff_ymid);
set(gca,'Box','on');
ylabel('CFF rate mid (MPa/yr)');
xlabel('\sigma_h / \sigma_v (-)');
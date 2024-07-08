

proj_folder = '/home/weng/Works/Softwares/PANTHER/';
folder      = sprintf('%s',[proj_folder,'/examples/Mmax/stress_data/']);
filenames   = dir(fullfile(folder,'*.mat'));
n_faults    = length(filenames);

global Mus Mud Dc pow_p;
Mus   = 0.6  ;
Mud   = 0.1  ;
Dc    = 0.01 ;
pow_p = 0.15 ;

for i = 1 : n_faults
    filename = filenames(i).name;
    load(sprintf('%s',[folder,filename]));
    disp(filenames(i).name)

    %% Read basic parameters
    y     = single_fault.pillars{1}.y   ; % Grid locations; Assumed the same for all pillars
    dy    = single_fault.pillars{1}.dy  ; % Grid length
    num_y = length(y)                   ; % Grid number
    num_p = single_fault.n_pillars      ; % number of pillars
    E     = single_fault.pillars{1}.input_parameters.young.value ;   % Young modulus; Assumed the same for all pillars
    nu    = single_fault.pillars{1}.input_parameters.poisson.value ; % Poisson ratio; Assumed the same for all pillars
    mu    = E / 2.0 / (1+nu)                         ;  % Shear modulus in mode III, unit is Pa
    mu_p  = mu/(1 - nu)                              ;  % Shear modulus in mode II, unit is Pa

    %% Initialize output
    G25D      = zeros(num_p);    % Elastic energy release: 0.5*Dtau*D
    Gc        = zeros(num_p);    % Fracture energy
    M0_unit   = 0;
    W         = zeros(num_p,3);    % Estimated rupture width
    Dtau      = zeros(num_p,num_y);    % Potential stress drop
    slip      = zeros(num_p,num_y);
    rup_state = -9; % rup_state = -1 means negative average stress drop
    % rup_state = 0  means there is negative slip
    % rup_state = 1  means rupture extention

    h1 = figure(1); clf(h1);

    for j = 1 : num_p
        thick = single_fault.pillars{j}.input_parameters.thick ;
        throw = single_fault.pillars{j}.input_parameters.throw ;
        sigma_long = single_fault.pillar_results{j}.stress{1}.sne(:,2) * 1e6 ;
        tau_long   = single_fault.pillar_results{j}.stress{1}.tau(:,2) * 1e6 - sigma_long*Mud;  % Potential stress drop
%        [top_ind,bot_ind,slip(j,:),Dtau(j,:),G25D(j),Gc(j)] = rupture_expend(tau_long,sigma_long,mu,dy);
    hold on
    plot(tau_long,'r-');
    end
    fig_name = sprintf('%s',[proj_folder,'/examples/Mmax/figs/',i,".png"]);
    saveas(h1,fig_name)
    close(1)

end


%% Functions
function [top_ind,bot_ind,slip_zeropad,Dtau_zeropad,G25D,Gc] = rupture_expend(tau_long,sigma_long,mu,dy)
top_ind = int16(length(tau_long)/2)-1;
bot_ind = int16(length(tau_long)/2)+1;

for i=1:min(length(tau_long)-bot_ind-1,top_ind-1)
    tau   = tau_long(top_ind:bot_ind);
    sigma = sigma_long(top_ind:bot_ind);
    tau_f = tau2tauf(tau,sigma,mu,dy);
    Dtau  = tau - tau_f;
    slip  = tau2slip(Dtau,dy,mu)  ;    
    Gc    = slip2Gc(slip, sigma)  ;
    [K_left,K_right] = tau2KII(Dtau,dy);
    G_left  = K_left^2/2/mu  ;
    G_right = K_right^2/2/mu ;

    if(G_left>Gc && G_right>Gc)
        top_ind = top_ind-1;
        bot_ind = bot_ind+1;
    elseif(G_left>Gc && G_right<=Gc)
        top_ind = top_ind-1;
    elseif(G_left<=Gc && G_right>Gc)
        bot_ind = bot_ind+1;
    else
        G25D = mean(0.5*Dtau.*slip) ;
        break
    end
end

slip_zeropad = zeros(length(tau_long),1);
Dtau_zeropad = zeros(length(tau_long),1);
slip_zeropad(top_ind:bot_ind) = slip(:);
Dtau_zeropad(top_ind:bot_ind) = Dtau(:);
end

%% Solve for final stress in 2D
function tau_f = tau2tauf(tau_0,sigma,mu,dy) 
global Mus Mud Dc pow_p;

num_p  = length(tau_0);
% Initialize D (starting point for iteration)
D = zeros(size(tau_0));  % Initial guess

% Set tolerance for convergence
tol = 1e-2            ;  % Example tolerance, adjust as needed
max_iterations = 1000 ;  % Maximum number of iterations (to prevent infinite loop)
total_pts = length(tau_0) ;
K   =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        K(n,m) = - mu/2/pi / dy / ((n-m)^2-1/4);
    end
end

for iter = 1:max_iterations

    F = tau_0 - K * D - f(D,sigma);  % find the zero
    J = - diag(f_prime(D,sigma)) - K ;  % Jacobian matrix

    % Solve for the update step using the given equation: tau_0 - K*D - f(D) 
    delta_D = J\F;
    D = D - 0.5* delta_D;

    for i=1:length(D)
        if(D(i)<0)
            D(i)=0;
        end
    end

    % Check convergence
    if norm(delta_D) < tol
        tau_f = f(D,sigma) ;
        return;
    end
end
tau_f = f(D,sigma) ;
disp(sprintf('The final norm of the iteration is %.2f',norm(delta_D)));
end


%% Caldualte the static stress intensity for asymmetric stress drop in 2D (Rice, 1980)
function [K_left,K_right] = tau2KII(Dtau,dx) %% Estimate G_II
epsilon = 0.0001 ;
num    =  length(Dtau) ;
W      =  (num-1)*dx ;
W_half =  W/2.0               ;
x      =  linspace(-W_half+epsilon,W_half-epsilon,num);
K_left =  sum((W_half-x).^0.5/(W_half+x).^0.5 .* Dtau) / (pi*W_half)^0.5;
K_right=  sum((W_half+x).^0.5/(W_half-x).^0.5 .* Dtau) / (pi*W_half)^0.5;
if(K_left<=0.0), K_left=0; end
if(K_right<=0.0), K_right=0; end
end

% Slip-stress relation in 2D
function slip = tau2slip(Dtau,dx,mu)
total_pts =  length(Dtau);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = - mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
stiff_inv = inv(stiff) ;
slip      = stiff_inv * Dtau;
end

function Dtau = slip2tau(slip,dx,mu)
total_pts =  length(slip);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = - mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
Dtau      = stiff *slip;
end

%%% Friction law and Fracture energy
function f = f(slip,sigma)
global Mus Mud Dc pow_p;
f=zeros(size(slip));
for i=1:length(slip)
    if(slip(i)<0)
        f(i) = sigma(i)*(Mus - Mud) ;
    else
        f(i) = sigma(i)*(Mus - Mud) / (1 + slip(i)/Dc/pow_p)^pow_p ;
    end
end
end

function f_prime = f_prime(slip,sigma)
global Mus Mud Dc pow_p;
f_prime=zeros(size(slip));
for i=1:length(slip)
    if(slip(i)<0)
        f_prime(i) = -sigma(i)*(Mus - Mud) / Dc ;
    else
        f_prime(i) = -sigma(i)*(Mus - Mud) / Dc / (1 + slip(i)/Dc/pow_p)^(pow_p+1) ;
    end
end
end

function Gc = slip2Gc(slip, sigma)
global Mus Mud Dc pow_p;
Gc=zeros(size(slip));
for i=1:length(slip)
    if(slip(i)<0)
        Gc = 0;
    else
        Gc = pow_p*Dc*sigma(i)*(Mus-Mud)*((1 + slip(i)/Dc)*(1+slip(i)/pow_p/Dc)^-pow_p-1)/(1-pow_p);
    end
end
end
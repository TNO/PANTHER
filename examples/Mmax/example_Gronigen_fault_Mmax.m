
proj_folder = '/home/weng/Works/Softwares/PANTHER/';
folder      = sprintf('%s',[proj_folder,'/examples/Mmax/stress_data/']);
filenames   = dir(fullfile(folder,'*.mat'));
n_faults    = length(filenames);

global Mus Mud Dc pow_p;
Mus   = 0.72  ;
Mud   = 0.22  ;
Dc    = 0.07 ;
pow_p = 4.0 ;

MM0  = cell(1,n_faults);
Mmax = cell(1,n_faults);

for i = 1 : 40
    filename = filenames(i).name;
    load(sprintf('%s',[folder,filename]));
    disp(filenames(i).name)

    %% Basic parameters and output variables
    y     = single_fault.pillars{1}.y   ; % Grid locations; Assumed the same for all pillars
    dy    = single_fault.pillars{1}.dy  ; % Grid length
    num_y = length(y)                   ; % Grid number
    num_p = single_fault.n_pillars      ; % number of pillars
    E     = single_fault.pillars{1}.input_parameters.young.value*1e6 ;   % Young modulus; Assumed the same for all pillars
    nu    = single_fault.pillars{1}.input_parameters.poisson.value ; % Poisson ratio; Assumed the same for all pillars
    mu    = E / 2.0 / (1+nu)                         ;  % Shear modulus in mode III, unit is Pa
    mu_p  = mu/(1 - nu)                              ;  % Shear modulus in mode II, unit is Pa

    G25D      = zeros(num_p,1);    % Elastic energy release: 0.5*Dtau*D
    Gc        = zeros(num_p,1);    % Fracture energy
    M0_unit   = 0;
    W         = zeros(num_p,3);    % Estimated rupture width
    Dtau      = zeros(num_p,num_y);    % Potential stress drop
    slip      = zeros(num_p,num_y);

    %% Check the potential shear stresses by plotting figures
    h1 = figure(1);
    subplot(2,1,1)
    for j = 1 : num_p
        sigma_long = single_fault.pillar_results{j}.stress{1}.sne(:,2) * 1e6 ;
        tau_long   = single_fault.pillar_results{j}.stress{1}.tau(:,2) * 1e6 - sigma_long*Mud;  % Potential stress drop
        hold on
        plot(tau_long,'r-');
    end
    subplot(2,1,2)
    hold on;
    plot(G25D,'k-');
    plot(Gc,'b-');
    fig_name = sprintf('%s',[proj_folder,'/examples/Mmax/figs/',filenames(i).name,".png"]);
    saveas(h1,fig_name)
    close(1)
% 
%     %% Solve for the rupture width for each pillar
%     for j = 1 : num_p
%         thick  = single_fault.pillars{j}.input_parameters.thick.value ;
%         throw  = single_fault.pillars{j}.input_parameters.throw.value ;
%         init_W = max(thick-throw,throw-thick) ;   % Assuming a minmum rupture width
%         sigma_long = single_fault.pillar_results{j}.stress{1}.sne(:,2) * 1e6 ;
%         tau_long   = single_fault.pillar_results{j}.stress{1}.tau(:,2) * 1e6 - sigma_long*Mud;  % Potential stress drop
%         if(max(tau_long)<0.0)  % Assumed as a strong barrier
%             G25D(j) = 1e-9;
%             Gc(j)   = 1 ;
%             continue ;         % No need to call the following function
%         end
%         [top_ind,bot_ind,slip(j,:),Dtau(j,:),G25D(j),Gc(j)] = rupture_width(tau_long,sigma_long,mu_p,dy,init_W);
%     end
% 
%     %% Solve for the Mmax for all pillars in one fault
%     MM0{i}  = calc_Mmax(single_fault,slip,G25D,Gc);
%     Mmax{i} = log10(MM0(:)*1e7)*2.0/3.0-10.7;
end



















%% Functions
%% Solve for M0 given an assumed nucleation pillar
function MM0 = calc_Mmax(single_fault,slip,G25D,Gc)

num_pillar  = length(G25D);
dy          = single_fault.pillars{1,1}.dy;
x_pillar    = zeros(1,num_pillar);
y_pillar    = zeros(1,num_pillar);
len_pillar  = zeros(1,num_pillar);
W_pillar    = zeros(1,num_pillar);
M0_pillar   = zeros(1,num_pillar);
Phi_pillar  = zeros(1,num_pillar);
GcG0_pillar = nan(1,num_pillar);
MM0         = nan(1,num_pillar);

for j = 1:num_pillar
    x_pillar(j) = [single_fault.pillar_info.Easting(j)];
    y_pillar(j) = [single_fault.pillar_info.Northing(j)];
    M0_pillar(j)  = sum(slip(j,:)*dy)*len_pillar(j);
    if(Gc(j)/G25D(j)>100.0 || Gc(j)/G25D(j)<0)
        Phi_pillar(j)  = NaN;
        W_pillar(j)    = NaN;
        GcG0_pillar(j) = NaN;
    else
        Phi_pillar(j) = (1-Gc(j)/G25D(j))*len_pillar(j);
        nonzeroslip   = find(slip(j,:)>0);
        W_pillar(j)   = (nonzeroslip(end)-nonzeroslip(1)+1)*dy;
        GcG0_pillar(j) = Gc(j)/G25D(j);
    end
end
% calculate along-strike length, M0, Phi of each pillarlen_fault
for j = 1:num_pillar-1
    len_pillar(j) = ((x_pillar(j) - x_pillar(j+1))^2 ...
        +(y_pillar(j) - y_pillar(j+1))^2 )^0.5;
end

for i=1:num_pillar
    if(isnan(GcG0_pillar(j)) || GcG0_pillar(j)>1)
        continue
    end
    [Rup_left, Rup_right] = rupture_length(i,len_pillar,W_pillar,GcG0_pillar);
    MM0(i) = sum(M0_pillar(Rup_left:Rup_right));
end
end

%% Solve for rupture length in horizontal profile (along stike)
function [Rup_left, Rup_right] = rupture_length(nuc_index,len,W,GcG0)
% Preallocate arrays
num   =  length(GcG0);
dx_W  = zeros(num);

for i=1:num
    if(W(i)>0.0)
        dx_W(i) = len(i)/W(i);
    else
        dx_W(i)=0.0;
    end
end
if (GcG0(nuc_index)>=1 || GcG0(nuc_index)<0) % cannot nucleate
    return
end
% Rupture to the left
for i = min(0,nuc_index-1):-1:1
    Rup_pot = sum((1-GcG0(i:nuc_index).*dx_W(i:nuc_index))) ;
    if (Rup_pot <= 0)
        Rup_left = i + 1;
        break;
    elseif(Rup_pot > 0 && i==1)
        Rup_left = 1;
        break;
    end
end

% Rupture to the right
for i = max(num,nuc_index+1):num
    Rup_pot = sum((1-GcG0(nuc_index:i).*dx_W(nuc_index:i))) ;
    if (Rup_pot <= 0)
        Rup_right = i - 1;
        break;
    elseif(Rup_pot > 0 && i==num)
        Rup_right = num;
        break;
    end
end
end


%% Solve for rupture width in 2D vertical profile (each pillar)
function [top_ind,bot_ind,slip_zeropad,Dtau_zeropad,G25D,Gc] = rupture_width(tau_long,sigma_long,mu,dy,init_W)
slip_zeropad = zeros(length(tau_long),1);
Dtau_zeropad = zeros(length(tau_long),1);
top_ind = int16(length(tau_long)/2)-int16(init_W/2);
bot_ind = int16(length(tau_long)/2)+int16(init_W/2);

for i=1:min(length(tau_long)-bot_ind-1,top_ind-1)
    tau   = tau_long(top_ind:bot_ind);
    sigma = sigma_long(top_ind:bot_ind);
    tau_f = tau2tauf(tau,sigma,mu,dy);     % Use Newton's method to solve the final friction strength
    if(max(tau_f)==0)
        G25D = 1e-9;  % Assumed as a strong barrier
        Gc   = 1 ;    % Assumed as a strong barrier
        return
    end
    Dtau  = tau - tau_f;
    [K_left,K_right] = tau2KII(Dtau,dy);
    slip  = tau2slip(Dtau,dy,mu)  ;        % Estimate slip based on the stress drop
    Gamma = slip2Gc(max(slip),max(sigma)); % Assuming the values of Gamma near top and bottom boundaries is controlled by the max slip
    G_left  = K_left^2/2/mu  ;
    G_right = K_right^2/2/mu ;
    if(G_left>Gamma && G_right>Gamma)
        disp("Double expension.")
        top_ind = top_ind-1;
        bot_ind = bot_ind+1;
    elseif(G_left>Gamma && G_right<=Gamma)
        disp("Left expension")
        top_ind = top_ind-1;
    elseif(G_left<=Gamma && G_right>Gamma)
        disp("Right expension")
        bot_ind = bot_ind+1;
    else
        Gc    = slip2Gc(slip, sigma)  ;
        G25D = mean(0.5*Dtau.*slip) ;
        break
    end
end

slip_zeropad(top_ind:bot_ind) = slip(:);
Dtau_zeropad(top_ind:bot_ind) = Dtau(:);
end

%% Caldualte the static stress intensity for asymmetric stress drop in 2D (Rice, 1980)
function [K_left,K_right] = tau2KII(Dtau,dx) %% Estimate G_II
epsilon = dx/2.0 ;
num    =  length(Dtau) ;
W      =  (num-1)*dx ;
W_half =  W/2.0             ;
x      =  zeros(size(Dtau)) ;
x(:)   =  linspace(-W_half+epsilon,W_half-epsilon,num);
K_left =  sum((W_half-x).^0.5./(W_half+x).^0.5.*Dtau*dx) /(pi*W_half)^0.5; % Caution! The format of Dtau and x shall be the same
K_right=  sum((W_half+x).^0.5./(W_half-x).^0.5.*Dtau*dx) /(pi*W_half)^0.5; % Caution! The format of Dtau and x shall be the same
if(K_left<=0.0), K_left=0; end
if(K_right<=0.0), K_right=0; end
end

%% Slip-stress relation in 2D
function slip = tau2slip(Dtau,dx,mu)
total_pts =  length(Dtau);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = - mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
slip      = stiff \ Dtau;
end

function Dtau = slip2tau(slip,dx,mu)
total_pts =  length(slip);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = - mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
Dtau      = stiff * slip;
end

%% Solve for final stress in 2D
function tau_f = tau2tauf(tau_0,sigma,mu,dy)
global Mus Mud Dc pow_p;

% Initialize D (starting point for iteration)
D = zeros(size(tau_0));
num_guess = 10000;
u = linspace(0,10,num_guess);
Res = zeros(num_guess,1);
for i=1:num_guess
    Res(i) = mean(tau_0) - (2*pi/4*mu/(dy*length(tau_0)))*u(i) - f(u(i),sigma);  % find the zero
end
[maxvalue,maxindex] = max(Res);
if(maxvalue<0)
    fprintf('This pillar is a strong barrier!\n');
    tau_f = 0;   % There is no solution for tau_f!
    return
else
    [~,minindex] = min(abs(Res(maxindex:end))) ;
    D(:) = u(minindex);   % An initial guess based on the average slip
end

% Set tolerance for convergence
tol = 1e-2            ;  % Example tolerance, adjust as needed
max_iterations = 10000 ;  % Maximum number of iterations (to prevent infinite loop)
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
    D = D -  delta_D;
    
    for i=1:length(D)
        if(D(i)<0)
            D(i)=0;
        end
    end

    % Check convergence
    if norm(delta_D) < tol
        tau_f = f(D,sigma) ;
        fprintf('The final norm is %.8f\n',norm(delta_D));
        return;
    end
end
tau_f = f(D,sigma) ;
fprintf('The final norm of the last iteration is %.8f\n',norm(delta_D));
end

%% Power-law friction law and Fracture energy
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

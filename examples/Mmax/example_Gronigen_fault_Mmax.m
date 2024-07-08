

proj_folder = '/home/weng/Works/Softwares/PANTHER/';
folder      = sprintf('%s',[proj_folder,'/examples/Mmax/stress_data/']);
filenames   = dir(fullfile(folder,'*.mat'));
n_faults    = length(filenames);




for i = 1 : 1
filename = filenames(i).name;
load(sprintf('%s',[folder,filename]));

%% Read basic parameters
y     = single_fault.pillars{1}.y   ; % Grid locations; Assumed the same for all pillars  
dy    = single_fault.pillars{1}.dy  ; % Grid length
num_y = length(y)                   ; % Grid number
num_p = single_fault.n_pillars      ; % number of pillars
E     = single_fault.pillars{1}.ensemble{1}.young ;   % Young modulus; Assumed the same for all pillars 
nu    = single_fault.pillars{1}.ensemble{1}.poisson ; % Poisson ratio; Assumed the same for all pillars 
mu    = E / 2.0 / (1+nu)                         ;  % Shear modulus in mode III, unit is Pa
mu_p  = mu/(1 - nu)                              ;  % Shear modulus in mode II, unit is Pa
eff_sigma = sigma(size(sigma)) 

%% Initialize output
Dtau      = zeros(num_p);    % Potential stress drop
G25D      = zeros(num_p);    % Elastic energy release: 0.5*Dtau*D
Gc        = zeros(num_p);    % Fracture energy
M0_unit   = 0;
W         = zeros(num_p,3);    % Estimated rupture width
slip      = zeros(num_p,num_y);
rup_state = -9; % rup_state = -1 means negative average stress drop
% rup_state = 0  means there is negative slip
% rup_state = 1  means rupture extention

for j = 1 : 1
    thick = single_fault.pillars{j}.ensemble{1}.thick ;
    throw = single_fault.pillars{j}.ensemble{1}.throw ;
    sigma = single_fault.pillar_results{j}.stress{1}.sne(:,1) ;
    tau_long = single_fault.pillar_results{j}.stress{1}.tau(:,1) ;
    [top_ind,bot_ind,slip(j,:),Dtau(j),G25D(j),Gc(j)] = rupture_expend(tau_long,eff_sigma,mu,dy);
end
end


% 
% %%%
% function [MM0,M0_seg,Rup_seg,Rup_ratio] = calc_Mmax(fault_data)
% 
% pillar_num = length(fault_data);
% x_fault    = zeros(1,pillar_num);
% y_fault    = zeros(1,pillar_num);
% 
% pillar_len = zeros(1,pillar_num);
% M0         = zeros(1,pillar_num);
% Phi        = zeros(1,pillar_num);
% 
% for j = 1:pillar_num
%     x_fault(j) = [fault_data{j}.x_coor];
%     y_fault(j) = [fault_data{j}.y_coor];
% end
% 
% % calculate along-strike length, M0, Phi of each pillar
% for j = 1:pillar_num-1
%     pillar_len(j) = ((x_fault(j) - x_fault(j+1))^2 ...
%         +(y_fault(j) - y_fault(j+1))^2 )^0.5;
%     M0(j)  = fault_data{j}.M0_unit.*pillar_len(j);
%     Phi(j) = (1-fault_data{j}.Gc/fault_data{j}.G0)*pillar_len(j);
% end
% Pillar_length = sum(pillar_len(1:pillar_num-1));
% 
% Allow_rupture = find(~isnan(Phi));
% [segment_num,segment_l,segment_r] = select_seg(Allow_rupture);
% M0_seg    = zeros(1,segment_num);
% Phi_seg   = zeros(1,segment_num);
% Rup_seg   = zeros(2,segment_num);
% Rup_ratio = zeros(1,segment_num);
% 
% for i = 1:segment_num
%     Phi_seg(i) = sum(Phi(segment_l(i):segment_r(i)));
%     if(Phi_seg(i)>=0)
%         M0_seg(i)    = sum(M0(segment_l(i):segment_r(i)));
%         Rup_seg(:,i) = [segment_l(i),segment_r(i)];
%         Rup_ratio(i) = (sum(pillar_len(segment_l(i):segment_r(i))))/Pillar_length;
%         continue
%     else
%         Positive_loc = find(Phi>=0);
%         Seg_pos = Positive_loc((Positive_loc>=segment_l(i) & Positive_loc<=segment_r(i)));
%         if(isempty(Seg_pos))
%             M0_seg(i)    =     0;
%             Rup_seg(:,i) = [0,0];
%             Rup_ratio(i) = 0;
%             continue
%         else
%             [nuc_num,nuc_l,nuc_r] = select_seg(Seg_pos);
%             M0_max   = 0.0;
%             for j=1:nuc_num
%                 for k=nuc_l(j):segment_r(i)
%                     Phi_nuc  = sum(Phi(nuc_l(j):k));
%                     if(Phi_nuc<0 || k==segment_r(i))
%                         M0_nuc  = sum(M0(nuc_l(j):k-1));
%                         break
%                     end
%                 end
%                 if(M0_nuc>M0_max)
%                     M0_seg(i)    =   M0_nuc;
%                     Rup_seg(:,i) =  [nuc_l(j),k-1];
%                     Rup_ratio(i) = (sum(pillar_len(Rup_seg(1,i):Rup_seg(2,i))))/Pillar_length;
%                     M0_max = M0_nuc;
%                 end
%             end
%             for j=1:nuc_num
%                 for k=nuc_r(j):-1:segment_l(i)
%                     Phi_nuc   = sum(Phi(k:nuc_r(j)));
%                     if(Phi_nuc<0 || k==segment_l(i))
%                         M0_nuc  = sum(M0(k+1:nuc_r(j)));
%                         break
%                     end
%                 end
%                 if(M0_nuc>M0_max)
%                     M0_seg(i)    =   M0_nuc;
%                     Rup_seg(:,i) =  [k,nuc_r(j)];
%                     Rup_ratio(i) = (pillar_len(Rup_seg(2,i))-pillar_len(Rup_seg(1,i)))/Pillar_length;
%                     M0_max = M0_nuc;
%                 end
%             end
%         end
%     end
% end
% MM0 = max(M0_seg);
% 
% %% Functions
%     function [segment_num,segment_l,segment_r] = select_seg(sequence)
%         segment_num   = 1;
%         segment_l(1)  = sequence(1);
%         for j = 2 : length(sequence)
%             if(sequence(j)-sequence(j-1)>1)
%                 segment_r(segment_num) = sequence(j-1);
%                 segment_num = segment_num + 1;
%                 segment_l(segment_num) = sequence(j);
%             end
%         end
%         segment_r(segment_num) = sequence(length(sequence));
%     end
% end


%% Functions
function [top_ind,bot_ind,slip,Dtau,G25D,Gc] = rupture_expend(tau_long,sigma,mu,dy)
top_ind = int16(length(tau_long)/2)-1;
bot_ind = int16(length(tau_long)/2)+1;

for i=1:min(length(tau_long)-bot_ind-1,top_ind-1)
    tau   = tau_long(top_ind:bot_ind);
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
end

%% Solve for final stress in 2D
function tau_f = tau2tauf(tau_0,sigma,mu,dy) 
% Initialize D (starting point for iteration)
D = zeros(size(tau_0));  % Initial guess

% Set tolerance for convergence
tol = 1e-6            ;  % Example tolerance, adjust as needed
max_iterations = 1000 ;  % Maximum number of iterations (to prevent infinite loop)
iter = 0;
total_pts = length(tau_0) ;
K   =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        K(n,m) = mu/2/pi / dy / ((n-m)^2-1/4);
    end
end
K_inv = inv(K) ;
while iter < max_iterations
    % Compute f(x)
    f_D = f(D,sigma) ;

    % Update x using the given equation: tau_0 - f(D) = K*D
    D_new = tau_0 - f_D;
    D_new = dot(K_inv,D_new); % Solve for x using matrix division (\)

    % Check convergence
    if norm(D_new - D) < tol
        D = D_new;
        break;
    end

    % Update x for the next iteration
    D = D_new;
    iter = iter + 1;
end
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

%% Slip-stress relation in 2D
function slip = tau2slip(Dtau,dx,mu)
total_pts =  length(Dtau);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
stiff_inv = inv(stiff) ;
slip      = dot(stiff_inv,Dtau);
end

function Dtau = slip2tau(slip,dx,mu)
total_pts =  length(slip);
stiff     =  zeros(total_pts,total_pts);
for n = 1:total_pts
    for m=1:total_pts
        stiff(n,m) = mu/2/pi / dx / ((n-m)^2-1/4);
    end
end
Dtau      = dot(stiff,slip);
end

%%% Friction law and Fracture energy
function f = f(slip,sigma)
Mus   = 0.6  ;
Mud   = 0.4  ;
Dc    = 0.01 ;
pow_p = 0.15 ;
size(slip)
size(sigma)
f = sigma.*(Mus - Mud) / (1 + slip./Dc./pow_p).^pow_p ;
end

function Gc = slip2Gc(slip, sigma)
Mus   = 0.6  ;
Mud   = 0.4  ;
Dc    = 0.01 ;
pow_p = 0.15 ;
Gc = pow_p*Dc*sigma*(Mus-Mud)*((1 + slip/Dc)*(1+slip/pow_p/Dc)^-pow_p-1)/(1-pow_p);
end


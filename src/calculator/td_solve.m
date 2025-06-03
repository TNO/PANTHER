% Pablo Ampuero
% Compute slip on "it" cells and stress on remaining cells
% induced by prescribed stress "itval" on "it" cells
% and slip "idval" on remaining cells.


% initial call on td_sovlve
% td_solve(icreep,0,-1,Ko) indices of creeping cells, itval; prescribed
% stress initially 0, 

% it: index of location where slip occurs, starting at location determined
% by min (tau
% idval   

% itval = excess shear stress in the slipping zone

function [t,d] = td_solve(it,itval,idval,K)

if length(itval) == 1, itval = repmat(itval,length(it),1); end
if length(itval) ~= length(it), error('it and itval must have same length'); end

[n,n2] = size(K);       % stiffness matrix
if n ~= n2, error('K must be a square matrix'); end

% location (index) of constrained-slip cells (non-"it" cells)
id = [1:n];     % all indices
id(it) = [];    % remove index of slipping cells it
if length(idval)==1, idval=repmat(idval,length(id),1); end
if length(idval)~=length(id), error('it and itval must have same length'); end

% Rearrange [stress]=-K*[slip] in an algebraic system A*x = y
% where x = unknown stresses and unknown slips,
% y = prescribed stresses and stress induced by prescribed slips
A = zeros(n,n);
A(:,it) = K(:,it);
At = zeros(n,1);
At(id) = 1;
A = A + diag(At);

y = zeros(n,1);
y = -K(:,id)*idval;     % non-slipping cells stress
y(it) = y(it) -itval;   % slipping cell stress

% Solve for x.
x = A\y;

% rearrange for output
t = zeros(n,1);
t(it) = itval;
t(id) = x(id);  % stress outside slipping cells
d = zeros(n,1);
d(id) = idval;
d(it) = x(it);  % displacement slipping cells

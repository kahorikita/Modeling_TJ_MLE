function y = get_trajectory(X,numofsim,Tmax,H,plant)
% simulate trajectory for feedback control model
%
% inputs:
%   X - parameters:
%       1. time of target jump
%       2. 
%       3. 
%       4-7. Feedback gains L
%   randomization seed
%   len ?
%   H ?
%   Ad, Bd - discretized plant model dynamics
Ad = plant.Ad;
Bd = plant.Bd;
delt = plant.delt;

rng(numofsim); % set seed for phase randomization

len = ceil(Tmax/delt); % total length of simulation
T = ceil(X(1)/delt); % timestep at which target jumps

L = X(4:7); % feedback gain matrix

% adjust high level parameters
% delt = 1/130; % time step length in secs

% parameters for A and B matrices
% G = params(1); I = params(2); tau = params(3);

% target = [-0.0001*ones(1,30) zeros(1,T(1)) -0.04*ones(1,T(2)-T(1)) 0.2*ones(1,180-T(2))];
target = [-0.0001*ones(1,50) zeros(1,T(1)) H*ones(1,len*2-T(1))]; %0.25
nstep = size(target,2);

% A = [0 1 0;0 -G/I 1/I;0 0 -1/tau]; % system dynamics matrix
% A2 = expm(delt*A); % discretize A
% Ad = blkdiag(A2,1); % augment matrix to include target 

% B = [0;0;1/tau;0]; % input matrix
% Bd = delt*B; % discretize B

order = size(Ad,1); % determine the order of the system

% initialize the state vector, x
x = zeros(order,nstep); % starting velocity and acceleration set to 0
x(1,1) = 0; % set starting hand position
x(4,1) = target(1,1); % set ending hand position

BdL = Bd*L;
% simulate trajectory
for i = 2:nstep
    x(:,i) = Ad*x(:,i-1) - BdL*x(:,i-1);
    % set target location
    x(4,i) = target(1,i);
end

% % simulate trajectory
% u = zeros(size(Bd,2),nstep); % movement commands
% sigma = 0.08; %0.05
% for i = 2:nstep
%     u(:,i) = -L*x(:,i-1);
% %     x(:,i) = Ad*x(:,i-1) + Bd*u(:,i) + Bd*sigma*randn(1); %add noise to u
%     x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);    
%     % set target location
%     x(4,i) = target(1,i);
% end

y.x = x(:,51:len+50);
y.T = T;
y.L = L;

% sampling = 130;
% figure(3); %cla
% hold on;
% xax2 = 1:length(y.x);
% time2 = 7.7*xax2;
% plot(time2,y.x(2,:),'b')
% load allsub_meanvel.mat;
% plot(mean_vel{2}(3,1:180)/(1/sampling)-mean_vel{2}(4,1:180)/(1/sampling)); % target on X, jump +y, MR
% axis([0 200 0 0.2]);





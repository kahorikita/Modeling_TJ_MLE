function sim = sim_vel_X1D_BSL(X,H,plant,Tmax)
Ad = plant.Ad;
Bd = plant.Bd;
delt = plant.delt;

Tjump = ceil(X(1)/delt); % mean time of Target Jump is gittered (gausian distribition)
sigma = abs(X(2))/delt; % variance of response time
% L = X(4:7); %one target jump
L = X(3:6); % feedback gain matrix

% adjust high level parameters
% delt = 1/130; % time step length in secs
start = 0; % starting position of the hand

% target = [zeros(1,T(1)) -0.04*ones(1,T(2)-T(1)) 0.2*ones(1,180-T(2))];
nstep = ceil(Tmax/delt)+ceil(sigma)*3; % number of timesteps to simulate
target = H*ones(1,nstep);
target(1:Tjump) = 0;

order = size(Ad,1); % determine the order of the system

% initialize the state vector, x
x = zeros(order,nstep); % starting velocity and acceleration set to 0
x(1,1) = start; % set starting hand position
x(4,1) = target(1,1); % set initial target position
% u = zeros(size(Bd,2),nstep); % movement commands

BdL = Bd*L;
Ad_BdL = Ad-BdL;
% simulate trajectory
for i = 2:nstep
    %u(:,i) = -L*x(:,i-1);
    %x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
    x(:,i) = Ad_BdL*x(:,i-1);
    % set target location
    x(4,i) = target(i);
end
% y_norm = normpdf(-ceil(X(2)*100):ceil(X(2)*100),0,X(3)*100); % gausian distribution

% convolve with Gaussian to account for variable response onset across
% trials
t_min = ceil(sigma)*3; % set convolution range to 3 standard deviations
y_norm = normpdf(-t_min:t_min,0,sigma); % gausian distribution

sim.x = x;
sim.convo = conv(x(2,:),y_norm,'same');
%sim.convo(:,len+1:end) = 0;

sim.T = Tjump;
sim.delt = delt;
sim.plant = plant;
sim.X = X;
sim.L = L;

% figure(tr); cla
% hold on;
% % xax1 = 1:length(sim.x(2,:));
% xax1 = 1:length(sim.convo);
% time1 = 7.7*xax1;
% % plot(time1,sim.x(2,:),'g')
% plot(time1,sim.convo,'g')
% xax2 = 1:length(y);
% time2 = 7.7*xax2;
% plot(time2,y,'r')
% 
% plot([7.7*T(1) 7.7*T(1)],[-0.1 0.4],'k');
% plot([7.7*T(2) 7.7*T(2)],[-0.1 0.4],'b');

% axis([0 1400 -0.05 0.25]);





function sim = sim_vel_X1D_MR(X,H,plant,Tmax)

Ad = plant.Ad;
Bd = plant.Bd;
delt = plant.delt;

Tjump = ceil(X(1)/delt); % mean time of Target Jump is gittered (gausian distribition) X(1):sec
Tjump2 = ceil(X(2)/delt);
sigma = abs(X(3))/delt; % variance of response time

start = 0; % starting position of the hand
nstep = ceil(Tmax/delt)+ceil(sigma)*3; % number of timesteps to simulate
target = -H*ones(1,nstep);
target(1:Tjump) = 0;
target(Tjump2:end) = H;

order = size(Ad,1); % determine the order of the system

% initialize the state vector, x
x = zeros(order,nstep); % starting velocity and acceleration set to 0
x(1,1) = start; % set starting hand position
x(4,1) = target(1,1); % set initial target position

R = X(4);
Q = [1 0 0 -1
    0 X(5) 0 0
    0 0 X(6) 0
    -1 0 0 1];
        
% calculate control law, L
n = 3000; % number of times to iterate over P 
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n
    P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
end
L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad); % control law


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

% convolve with Gaussian to account for variable response onset across
% trials
t_min = ceil(sigma)*3; % set convolution range to 3 standard deviations
y_norm = normpdf(-t_min:t_min,0,sigma); % gausian distribution

sim.x = x;
sim.convo = conv(x(2,:),y_norm,'same');
sim.acc = diff(sim.convo)/delt;
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





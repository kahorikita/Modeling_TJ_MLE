
%% target not jump/BSL/MR 
clear all;
%load TJ_BSL_nojump_targ_on_x;
%load TJ_BSL_nojump_targ_on_y;
load TJ_MR_nojump_targ_on_x;
load TJ_MR_nojump_targ_on_y;
time = 1:260;
time = time*7.7;
Hz = 130.004;
rng(2);

% adjust high level parameters
delt = 0.001; % time step length in secs
start = [0 0]; % starting position of the hand

% Single joint reaching movements:
G = 0.14;        % Viscous Constant: Ns/m; originally 0.14
I = 0.1;         % Inertia Kgm2; originally 0.1
tau = 0.066;    % Muscle time constant, s; originally 0.066

for targ=1:2
    
    if targ==1 % target on y axis
        target = [zeros(1,1500)
            0.08*ones(1,1500)];
    elseif targ==2 % target on x axix
        target = [0.08*ones(1,1500)
            zeros(1,1500)];
    end
    
    nstep = size(target,2);
    
    % create state space model in discrete time
    A = [0 1 0
        0 -G/I 1/I
        0 0 -1/tau]; % system dynamics matrix
    B = [0 0 1/tau]';
    
    z = size(A,1)+1;
    Ad = expm(A*delt);
    Ad = blkdiag(Ad, 1);
    Ad = [Ad zeros(z)
        zeros(z) Ad];
    
    Bd = delt*B;
    Bd = [Bd;0];
    Bd = [Bd zeros(z,1)
        zeros(z,1) Bd];
    %BSL
%     R = [0.002 0
%         0 0.002]; % effort cost- default is 0.0001
    %MR
    R = [0.002 0
        0 0.002];
    
    Q2 = [1 0 0 -1
        0 0.04 0 0
        0 0 0.01 0
        -1 0 0 1];
    Q = [Q2 zeros(z)
        zeros(z) Q2];
    
    % state vector
    order = size(Ad,1); % order of the system
    x = zeros(order,nstep);
    x(1,1) = start(1); % x position
    x(end/2+1,1) = start(2); % y position
    x(end/2,1) = target(1,1); % starting target x position
    x(end,1) = target(2,1); % starting target y position
    u = zeros(size(Bd,2),nstep); % movement commands
    
    % calculate feedback gain
    n = 5000;
    P = zeros(order,order,n);
    P(:,:,1) = rand(order);
    for i = 2:n
        P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
    end
    L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad)*.2;
    
    % simulate trajectory
    for i = 2:nstep
        u(:,i) = -L*x(:,i-1);
        x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
        
        % set target location
        x(end/2,i) = target(1,i);
        x(end,i) = target(2,i);
    end
    
    figure(1); 
    subplot(1,2,targ)
    hold on;
    plot(x(1,:),x(5,:)) % hand trajectory
    plot(x(4,1),x(8,1), 'o') % target starting position
    plot(x(4,end),x(8,end), 'o') % target end position
    axis([-0.1 0.1 0 0.1])
    legend('Hand trajectory','Target Start','Target End','Location','northwest')
    
    figure(2);
    subplot(1,2,targ)
    hold on;
    if targ==1
        plot(x(6,:),'LineWidth',1.5) % when target on y
        title('target on y, y veloicty')
        plot(time, mean_vel_y,'LineWidth',1.5)
        ylabel('Velocity (m/s)')
        xlabel('Time (s)')
    elseif targ==2
        plot(x(2,:),'LineWidth',1.5) % when target on y
        title('target on x, x veloicty')
        plot(time, mean_vel_x,'LineWidth',1.5)
        ylabel('Velocity (m/s)')
        xlabel('Time (s)')
    end
    
end


%% target jump / BSL
clear all;
%load dat.mat 
%load angle.mat
load allsub_meanvel.mat
time = 1:259;
time = time*7.7;
sampling = 130.004;
rng(2);

% adjust high level parameters
delt = 0.001; % time step length in secs
start = [0 0]; % starting position of the hand

% Single joint reaching movements:
G = 0.14;        % Viscous Constant: Ns/m; originally 0.14
I = 0.1;         % Inertia Kgm2; originally 0.1
tau = 0.066;    % Muscle time constant, s; originally 0.066

% target trajectory; row 1: x position, row 2: y position
% target = [zeros(1,1000) ones(1,2000)
%           ones(1,3000)];
% target = [zeros(1,300) 0.05*ones(1,1200)
%             0.08*ones(1,1500)];

% 1: target on Y axsis, target jump +x
% 2: target on Y axsis, target jump -x
% 3: target on X axsis, target jump +y
% 4: target on X axsis, target jump -y

for targ = 1:4
    
    if targ == 1
        % target on Y, jump +x
        target = [zeros(1,260) 0.05*ones(1,1240)
            0.08*ones(1,1500)];
    elseif targ == 2
        % target on Y, jump -x
        target = [zeros(1,260) -0.05*ones(1,1240)
            0.08*ones(1,1500)];
    elseif targ == 3
        % target on X, jump +y
        target = [0.08*ones(1,1500)
            zeros(1,260) 0.05*ones(1,1240)];
    elseif targ == 4
        % target on X, jump -y
        target = [0.08*ones(1,1500)
            zeros(1,260) -0.05*ones(1,1240)];
    end
    nstep = size(target,2);
    
    % create state space model in discrete time
    A = [0 1 0
        0 -G/I 1/I
        0 0 -1/tau]; % system dynamics matrix
    B = [0 0 1/tau]';

    z = size(A,1)+1;
    Ad = expm(A*delt);
    Ad = blkdiag(Ad, 1);
    Ad = [Ad zeros(z)
        zeros(z) Ad];
    
    Bd = delt*B;
    Bd = [Bd;0];
    Bd = [Bd zeros(z,1)
        zeros(z,1) Bd];
    
    % accuracy and effort costs
    % R = [0.0001 0
    %     0 0.0001]; % effort cost- default is 0.0001
    % Q2 = [1 0 0 -1
    %     0 0.1 0 0
    %     0 0 0 0
    %     -1 0 0 1];
    %%% 04/23/19 probably better to change only effort costs(R) and
    %%% accerelation
    %%% best params for target jump
%     R = [0.0035 0
%         0 0.0035]; % effort cost- default is 0.0001
%     Q2 = [1 0 0 -1
%         0 0.01 0 0
%         0 0 0.005 0
%         -1 0 0 1];
%     Q = [Q2 zeros(z)
%         zeros(z) Q2];
    
    R = [0.002 0
        0 0.002]; % effort cost- default is 0.0001
    Q2 = [1 0 0 -1
        0 0.04 0 0
        0 0 0.01 0
        -1 0 0 1];
    Q = [Q2 zeros(z)
        zeros(z) Q2];
    
    % state vector
    order = size(Ad,1); % order of the system
    x = zeros(order,nstep);
    x(1,1) = start(1); % x position
    x(end/2+1,1) = start(2); % y position
    x(end/2,1) = target(1,1); % starting target x position
    x(end,1) = target(2,1); % starting target y position
    u = zeros(size(Bd,2),nstep); % movement commands
    
    % calculate feedback gain
    n = 5000;
    P = zeros(order,order,n);
    P(:,:,1) = rand(order);
    for i = 2:n
        P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
    end
    L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad)*.2;
    
    % simulate trajectory
    for i = 2:nstep
        u(:,i) = -L*x(:,i-1);
        x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
        
        % set target location
        x(end/2,i) = target(1,i);
        x(end,i) = target(2,i);
    end

    figure(1); hold on;
    %subplot(1,4,targ); hold on;
    plot(x(1,:),x(5,:)) % hand trajectory
    plot(x(4,1),x(8,1), 'o') % target starting position
    plot(x(4,end),x(8,end), 'o') % target end position
    %axis([-1 1 0 1])
    legend('Hand trajectory','Target Start','Target End','Location','northwest')
    
    
    figure(2);
    subplot(1,4,targ); hold on;
    if targ == 1
        plot(x(2,:),'LineWidth',1.5) % when target on y
        title('target on y, jump +x')
    elseif targ == 2
        plot(x(2,:),'LineWidth',1.5) % when target on y
        title('target on y, jump -x')
    elseif targ == 3
        plot(x(6,:),'LineWidth',1.5) % when targe on x
        title('target on x, jump +y')
    elseif targ == 4
        plot(x(6,:),'LineWidth',1.5) % when targe on x
        title('target on x, jump -y')
    end
    plot(time, mean_vel{1}(targ,:)/(1/sampling),'LineWidth',1.5)
    ylabel('Velocity (m/s)')
    xlabel('Time (s)')
        
end





%% target jump / mirror

clear all;
load allsub_meanvel.mat
time = 1:259;
time = time*7.7;
sampling = 130.004;
rng(2);

% adjust high level parameters
delt = 0.001; % time step length in secs
start = [0 0]; % starting position of the hand

% Single joint reaching movements:
G = 0.14;        % Viscous Constant: Ns/m; originally 0.14
I = 0.1;         % Inertia Kgm2; originally 0.1
tau = 0.066;    % Muscle time constant, s; originally 0.066

% 1: target on Y axsis, target jump +x
% 2: target on Y axsis, target jump -x
% 3: target on X axsis, target jump +y
% 4: target on X axsis, target jump -y
% target trajectory; row 1: x position, row 2: y position

% target on X --> (jump -y) --> jump +y
target = [0.08*ones(1,1700)
    zeros(1,260) -0.025*ones(1,430-260) 0.05*ones(1,1700-(430-260)-260)];
nstep = size(target,2);

% create state space model in discrete time
A = [0 1 0
    0 -G/I 1/I
    0 0 -1/tau]; % system dynamics matrix
B = [0 0 1/tau]';

z = size(A,1)+1;
Ad = expm(A*delt);
Ad = blkdiag(Ad, 1);
Ad = [Ad zeros(z)
    zeros(z) Ad];

Bd = delt*B;
Bd = [Bd;0];
Bd = [Bd zeros(z,1)
    zeros(z,1) Bd];

% accuracy and effort costs
% R = [0.0001 0
%     0 0.0001]; % effort cost- default is 0.0001
% Q2 = [1 0 0 -1
%     0 0.1 0 0
%     0 0 0 0
%     -1 0 0 1];
%%% 04/23/19 probably better to change only effort costs(R) and
%%% accerelation
%%% best params for mirror
% R = [0.0035 0
%     0 0.0035]; % effort cost- default is 0.0001
% Q2 = [1 0 0 -1
%     0 0.1 0 0
%     0 0 0.02 0
%     -1 0 0 1];
% Q = [Q2 zeros(z)
%     zeros(z) Q2];

R = [0.002 0
    0 0.002]; % effort cost- default is 0.0001
Q2 = [1 0 0 -1
    0 0.04 0 0
    0 0 0.01 0
    -1 0 0 1];
Q = [Q2 zeros(z)
    zeros(z) Q2];


% state vector
order = size(Ad,1); % order of the system
x = zeros(order,nstep);
x(1,1) = start(1); % x position
x(end/2+1,1) = start(2); % y position
x(end/2,1) = target(1,1); % starting target x position
x(end,1) = target(2,1); % starting target y position
u = zeros(size(Bd,2),nstep); % movement commands

% calculate feedback gain
n = 5000;
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n
    P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
end
L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad)*.2;

% simulate trajectory
for i = 2:nstep
    u(:,i) = -L*x(:,i-1);
    x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
    
    % set target location
    x(end/2,i) = target(1,i);
    x(end,i) = target(2,i);
end

figure(1); clf; hold on;
plot(x(1,:),x(5,:)) % hand trajectory
plot(x(4,1),x(8,1), 'o') % target starting position
plot(x(4,end),x(8,end), 'o') % target end position
%axis([-1 1 0 1])
legend('Hand trajectory','Target Start','Target End','Location','northwest')

figure(2); clf; hold on;
plot(x(6,:),'LineWidth',1.5) % when targe on x
plot(time, mean_vel{2}(3,:)/(1/sampling),'LineWidth',1.5)
title('target on x, jump +y')
ylabel('Velocity (m/s)')
xlabel('Time (s)')




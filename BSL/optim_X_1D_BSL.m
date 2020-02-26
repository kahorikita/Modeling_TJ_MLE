%%% Linit is calcuted by optim_L_1D

%%% -define the true parameters
%%% -simulate model to get 'true' trajectory
%%% -save trajectory to compare with candidate model fits in get_error function
clear
clc

tic

cols(:,:,1) = [ 0 210 255; 255 210 0; 0 0 0; 210 0 255]/256;
cols(:,:,2) = [ 0 155 255; 255 100 0; 0 0 0; 155 0 255]/256;
cols(:,:,3) = [ 0 100 255; 255 0 0; 0 0 0; 100 0 255]/256;

% for simulate velocity 
params = [0.14 0.1 0.066]; % Single joint reaching movements:
G = params(1); % Viscous Constant: Ns/m; originally 0.14    ------%%% - Can you provide a ref for where these values come from?
I = params(2); % Inertia Kgm2; originally 0.1
tau = params(3); % Muscle time constant, s; originally 0.066
delt = .001; % time step length in secs
A = [0 1 0;0 -G/I 1/I;0 0 -1/tau]; % system dynamics matrix
A = blkdiag(A,0); % add target state
B = [0;0;1/tau;0]; % input matrix

% discretize (use matrix exponential trick to compute Ad and Bd in single
% timestep (see 'Discretization' on Wikipedia)
M = expm([A B;zeros(1,5)]*delt);
Ad = M(1:4,1:4);
Bd = M(1:4,5);
%TrueParam = [35/100 50/100 1.56 0.6 -0.06 -0.87];

plant.A = A;
plant.B = B;
plant.delt = delt;
plant.Ad = Ad;
plant.Bd = Bd;

Tmax = .25; % max. time of simulation in s
len = ceil(Tmax/delt);
numofbootstraps = 1; % number of bootstraps
H = 0.1; % amplitude of target, displacement: -0.05, +0.05 -->0.1

% params for initial value of target jump
istart = 30; %20
numofsteps = 1; %15
Hz = 130.004;

%%

% initilization
Result_all{numofbootstraps} = [];

% switch data used in simulation
% simulate velocity / mean of all subs / bootstrap for each sub
disp('1:simulated velocity, 2:mean of all subs, 3:bootstrap data, 4:mean of each sub');  
m = input('Choose data: ');
switch m
    case 1 % to simulate velocity profiles for parameter recovery
        disp('simulated velocity')  
%         TrueParam = Result{1,1};
        TrueParams = [.15 .025 1.56 0.6 0.06 -0.87]; %[TJ, 2 params of gausian distribution, 4 params of L] 
        %for numofsim = 1:100 % it's different from the "tr" below   
            %y_temp = get_trajectory(TrueParam,numofsim,Tmax,H,plant);
            sim = sim_vel_X1D_BSL(TrueParams,H,plant,Tmax)
            %vel(numofsim,:) = y_temp.x(2,:); 
        %end
        y = sim.convo;
        numofbootstraps = 1;
        numofsubjects = 1;
    case 2 % average velocity of all subjects (experiment data)
        disp('mean of all subs')
        load allsub_meanvel.mat
%         y = mean_vel{2}(3,1:len)*Hz-mean_vel{2}(4,1:len)*Hz; % target on X, jump +/-y, MR
        y_130 = mean_vel{1}(3,:)-mean_vel{1}(4,:); % target on X, jump +/-y, BSL
        % resample to 1000Hz
        y = resample(y_130,1000,130)*130;
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        
        numofbootstraps = 1;
        numofsubjects = 1;
    case 3 % each subject data with bootstrap
        disp('bootstrap data')
        load AllData.mat
        numofsubjects = 20;
    case 4 % mean of each subject data
        disp('mean of each sub')
        load AllData.mat
        numofsubjects = 20;
    case 5 % average acceleration of all subjects (experiment data)
        disp('mean acceleration of all subs')
        load allsub_meanvel.mat
        y_130 = diff(mean_vel{1, 1}(3,:)*130)*130-diff(mean_vel{1, 1}(4,:)*130)*130; % target on X, jump +/-y, BSL
        % resample to 1000Hz
        y = resample(y_130,1000,130);
        y = y(101:end); % subtract 100 ms instrument delay
        y = y - mean(y(1:100)); % subtract baseline
        
        numofbootstraps = 1;
        numofsubjects = 1;
    otherwise
        disp('other value')
end

disp('1:MLE, 2:BADS');
optm = input('Choose an optimization method: '); %
switch optm
    case 1
        disp('MLE')
    case 2
        disp('BADS')
end

c = 1;
fhandle = figure(c); clf; hold on
% set(fhandle, 'Position', [200*c, 100, 450, 450]); % set size and loction on screen
set(fhandle, 'Position', [200*c, 100, 900, 450]); % set size and loction on screen
set(fhandle, 'Color','w') % set background color to white 
set(gca,'FontSize',10);

for subject = 1:numofsubjects
% for subject = 2:2
    
    if m == 3
        datap = data(subject,3+2).Vel_CrX_post*Hz; % target jump +y data{1,3} data(subject,1:2)=TR
        datam = data(subject,4+2).Vel_CrX_post*Hz; % target jump -y data{1,4}
        ybs = get_bootstrap(datap,datam,len,numofbootstraps);
    elseif m == 4
        datap = data(subject,3+2).Vel_CrX_post*Hz;
        datam = data(subject,4+2).Vel_CrX_post*Hz;
        ytmp = nanmean(datap) - nanmean(datam);
        stable = nanmean(ytmp(1:13));
        ytmp = ytmp-stable; % velocity starts from zero
        ytmp = ytmp(1:len);
        y = ytmp;   
    end
    
    for tr = 1:numofbootstraps % number of bootstrap
        
        Xopt = [];
        
        if m == 3
            % set one bootstrap data
            y = ybs(tr,:);
        end
        
        for i = 1%istart:istart+numofsteps
            
%             TrueParam = [i/100 0.04 0.15 1.56 0.6 -0.06 -0.87]; % target jump, µ, sigma, L
%             Tinit = TrueParam(1); %one target jump
%             Ginit = TrueParam(2:3);
%             Linit = TrueParam(4:7);
%             Xinit = [Tinit Ginit Linit];
%             target = [zeros(1,ceil(Tinit(1)*100)) 0.3*ones(1,len-ceil(Tinit(1)*100))];
            %TrueParam = [.1 .01 0.6 -0.06 -0.87]; % [jumpTime 
            
            Ginit = [.2 .1]; % initial guess for Gaussian parameters
            Linit = [1 0.5 0.1 -1] % initial guess for feedback gains
            
            %Ginit = TrueParams(1:2);
            %Linit = TrueParams(3:6);
            Xinit = [Ginit Linit];
            %target = [zeros(1,ceil(Ginit(1)*100)) 0.3*ones(1,len-ceil(Ginit(1)*100))];
            
            f_targ = @(X) get_error_X1D_BSL(X,y,H,plant,Tmax);
            
            
            switch optm
                case 1 % fmincon
                    Aeq = [];
                    beq = [];
                    % used for TA and myu
%                     lb = [0.25 0.001 0.001 0 0 -2 -2]; % Lower
%                     ub = [0.4 1 1 5 1 1 0]; % Upper bounds
                    lb = [0 0.001 -2 -2 -2 -2]; % Lower bounds
                    ub = [0.5 0.5 2 2 2 2]; % Upper bounds
                    
                    Xopt = fmincon(f_targ,Xinit,[],[],Aeq,beq,lb,ub);
                    %Xopt = fmincon(f_targ,Xinit);
                    %Xopt = fminsearch(f_targ,TrueParams)
                    
                case 2 % bads
%                     plb = [0.25 0.001 0.001 0 0 -2 -2]; % Lower bounds
%                     pub = [0.4 1 1 5 1 1 0]; % Upper bounds
%                     lb = [0 0.001 0.001 0 0 -2 -2]; % Plausible lower bounds
%                     ub = [0.4 1 1 5 1 1 0]; % Plausible upper bounds
                    lb =  [0   0.001 0  0  0 -5]; % Lower bounds                    
                    plb = [0.1 0.005 .7 0 0 -2]; % Plausible Lower bounds
                    
                    pub = [0.3 0.05  2  2  1  0]; % Plausible Upper bounds
                    ub =  [0.5 .2    5  2  1  0]; % Upper bounds
                    
                    Xopt = bads(f_targ,Xinit,lb,ub,plb,pub);
                otherwise
            end
            
            %ferror = get_error_X1D_BSL(temp,y,len,H,Ad,Bd);
            %Xopt = [Xopt; temp ferror];
            
        end
        %Result_all{subject,tr} = Xopt;
        
        % plot optimal result = minimum error
        %[M, I] = min(Result_all{subject,tr}(:,end)); %one target jump
        opt = sim_vel_X1D_BSL(Xopt,H,plant,Tmax);
        %--HACK--
        %if(size(opt.convo)>len)
        %    opt.convo = opt.convo(:,1:len);
        %end
        
        %Resultopt{subject}.param(tr,:) = Result_all{subject,tr}(I,:);
        %Resultopt{subject}.trajectory(tr,:) = opt.convo;
        %Resultopt{subject}.I(tr,:) = I;
        %Resultopt{subject}.mean = mean(Resultopt{subject}.param(:,1));
        %Resultopt{subject}.var = var(Resultopt{subject}.param(:,1));
        
%         figure(subject);
%         subplot(5,10,tr); hold on;
        
        figure(1); clf; hold on;
        %title(['participant ',num2str(subject)],'fontsize',10); %,model(m).name,' model'
        %x = 7.7*(1:len)-100; % subtract 100ms
        plot(y,'color',cols(4,:,c),'linewidth',1.5)
%         plot(x,opt.convo(:,1:len),'g','linewidth',1.1)
        plot(opt.convo(:,1:len),'color',cols(1,:,c),'linewidth',1.5)
        legend('simulated data','fit')
        %plot([7.7*100*Result_all{subject,tr}(I,1)-100 7.7*100*Result_all{subject,tr}(I,1)-100],[-0.1 0.4],'k','linewidth',1);
        xlabel('Time (ms)','FontSize',10)
        ylabel('Velocty (m/s)','FontSize',10)
        %axis([0 700 -0.1 0.5])
        
        
    end
end

% save 'TJModelFits.mat' Resultopt Result_all


toc

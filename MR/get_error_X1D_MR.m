function e2 = get_error_X1D_MR(X,y,H,Hm,plant,Tmax,m)

Tmax_sim = (size(y,2))*plant.delt + .5; % always simulate 500 ms further than the data

% perform simulations
sim = sim_vel_X1D_MR(X,H,Hm,plant,Tmax_sim); 

imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data
if(imax>size(y,2))
    error('Tmax too long')
end

length_sim = size(sim.x,2);
length_data = size(y,2);
%--HACK--
%if(size(sim.convo,2)>len) %%%180
%    sim.convo = sim.convo(:,1:len);
%end

% calculate error between simulation and acual data 
if m == 5 || m == 6 % input is acceleration
    e2 = nanmean((y(1:imax)-sim.acc(1:imax)).^2);  
    figure(3); clf; hold on
    plot(y,'m')
    plot(sim.acc,'b')
    plot((y(1:imax)-sim.acc(1:imax)),'r')

else % input is velocity
    e2 = nanmean((y(1:imax)-sim.x(1:imax)).^2);
    figure(3); clf; hold on
    plot(y,'m')
    plot(sim.x,'b')
    plot((y(1:imax)-sim.x(1:imax)),'r')
end


% % calculate error between simulation and acual data 
% % e2 = nanmean((y-sim.x(2,:)).^2);
% e2 = nanmean((y(1:imax)-sim.convo(1:imax)).^2);
% % [tr X]


% disp(num2str(e2))
%
figure(3); clf; hold on
plot(y,'m')
plot(sim.convo,'b')
plot((y(1:imax)-sim.convo(1:imax)),'r')
%}


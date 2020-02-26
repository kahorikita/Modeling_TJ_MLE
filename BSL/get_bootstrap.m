function y = get_bootstrap(datap,datam,len,numofbootstrap)

for i=1:numofbootstrap
    % 1. bootstrap and calculate mean velocity
    % for target jump +y
    % datap = data{1,3}.Vel_CrX_post*Hz;
    bootp = datasample(datap,size(datap,1));
    
    % 2. for target jump -y
    % datam = data{1,4}.Vel_CrX_post*Hz;
    bootm = datasample(datam,size(datap,1));
    
    % 3. calculate diff_vel of mean velocities
    ytmp = nanmean(bootp) - nanmean(bootm);
    stable = nanmean(ytmp(1:13));
    ytmp = ytmp-stable; % velocity starts from zero
    ytmp = ytmp(1:len);
    y(i,:) = ytmp;
end

% 4. optimization

% 5. repeat 1-4


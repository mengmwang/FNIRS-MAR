% TestExperimental - Test the motion artefact removal algorithm using
% experimental fNIRS signals and simulated artefacts


%% define the functions
algo = mar_algo;

%% load data
% experimental fNIRS signal 
load('experimental_signal.mat');

% x_data: one channel experimental fNIRS signal
% yma: experimental signal with simulated motion artefacts


%% define parameters
% parameters can be tuned for better results
W = 200;  % window size
inc = 100; % incremental window size
r = 0.2;
coeff_r = round(W*r);
al_w = 0.1;

%% algorithm
for i = 1:1:(size(yma,1)-W)/inc+1
    yma_buffer(:,i)=yma((i-1)*inc+1:(i-1)*inc+W);
end

if W == size(yma,1)
    yma_buffer=yma;
end

y_theta_w_buf = zeros(size(yma_buffer));

% for loop
parfor p = 1:1:size(yma_buffer,2)
    yma_t = yma_buffer(:,p);
    yt = algo.diff(yma_t,1);  % calculate differences

    [br,brt] = algo.rb(yma_t,coeff_r);  % reduced basis functions

    % robust estimation
    theta_ls = algo.ls(brt,yt); 
    theta_w = algo.estimate(brt,yt,theta_ls,al_w);
    y_theta_w_t = br*theta_w;

    y_theta_w_buf(:,p) = y_theta_w_t;
end

y_theta_w_all = zeros(size(yma,1),size(yma_buffer,2));

for j =1:1:size(yma_buffer,2)
    y_theta_w_all((j-1)*inc+1:(j-1)*inc+W)=y_theta_w_buf(:,j);
end

y_theta_w = sum(y_theta_w_all,2)./sum(y_theta_w_all~=0,2);

%% 
% plot 
figure(1)
subplot(211)
plot(x_data)
hold on
plot(yma)
legend('Experimental Signal','Signal with Simulated Motion Artefacts')
subplot(212)
plot(x_data)
hold on
plot(y_theta_w)
legend('Experimental Signal','After Artefact Removal Algorithm')
suptitle('Motion artefact removal on experimental signal (simulated artefacts)')

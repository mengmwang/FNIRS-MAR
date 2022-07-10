% TestSimulation - Test the motion artefact removal algorithm using
% simulated signals and simulated artefacts


%% define the functions
algo = mar_algo;

%% load data
% load simulated data

load('simulated_signal.mat');
N = size(x_data,1);
n = (1:N)';

% x_data: low frequency components + random noise
% yma: signal with simulated motion artefacts

%% define parameters and algorithm
% parameters can be tuned
coeff_r = round(N*.15);
al_w = 0.6;

% calculate the first order difference
yt = algo.diff(yma,1);
% form reduced basis functions
[br,brt] = algo.rb(yma,coeff_r);

% robust estimation
theta_ls = algo.ls(brt,yt);
theta_w = algo.estimate(brt,yt,theta_ls,al_w);
y_theta_w = br*theta_w;


%% 
% plot 
figure(1)
subplot(211)
plot(x_data)
hold on
plot(yma)
legend('Simulated Signal','Signal with Motion Artefacts')
subplot(212)
plot(x_data)
hold on
plot(y_theta_w)
legend('Simulated Signal','After Artefact Removal Algorithm')
suptitle('Motion artefact removal on simulated signal (simulated artefacts)')


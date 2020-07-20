clear;clc;  close all;
addpath(genpath('TestData'));

% Test Dataset IDs
ID = { 'DATA_01_TYPE01', 'DATA_02_TYPE02', 'DATA_03_TYPE02', 'DATA_04_TYPE02', ...
    'DATA_05_TYPE02', 'DATA_06_TYPE02', 'DATA_07_TYPE02', 'DATA_08_TYPE02',...
    'DATA_09_TYPE02', 'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02'};
IDt = { 'DATA_01_TYPE01_BPMtrace', 'DATA_02_TYPE02_BPMtrace', 'DATA_03_TYPE02_BPMtrace', 'DATA_04_TYPE02_BPMtrace', ...
    'DATA_05_TYPE02_BPMtrace', 'DATA_06_TYPE02_BPMtrace', 'DATA_07_TYPE02_BPMtrace', 'DATA_08_TYPE02_BPMtrace',...
    'DATA_09_TYPE02_BPMtrace', 'DATA_10_TYPE02_BPMtrace','DATA_11_TYPE02_BPMtrace','DATA_12_TYPE02_BPMtrace'};

% creating Adaptive Filter object
M = 32;                    % Filter order
lam = 0.9999;              % RLS forgetting factor
delta = 0.1;               % Initial input covariance estimate
w0 = zeros(M,1);           % Initial tap weight vector
P0 = (1/delta)*eye(M,M);   % Initial setting for the P matrix
Zi = zeros(M-1,1);         % FIR filter initial states
ha = adaptfilt.rls(M,lam,P0,w0,Zi);
f=linspace(0,4,300); % frequency range

for idnb = 1:12
    
    load(ID{idnb});                          % load test dataset
    load(IDt{idnb}); 
    
    srate = 125;                             % Fs=125 Hz
    
    window   = 8 * srate;                    % window length is 8 seconds
    step     = 2 * srate;                    % step size is 2 seconds
    
    windowNb = (length(sig)-window)/step + 1;  % total number of windows(estimates)
    
    %**********************************************************************
    % Please write your codes as follows (i.e.,inputing data window by window)
    for i =   1  :  windowNb
        
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        
        % Our algorithm's code
        
        dl=12; % No. of samples to be delayed in acceleromter data
        
        if i==1
            prev=0; % for the 1st window, no previous data are available
            bpm(1)=0;
            ty=10; % flag variable
            tf=2; % flag variable
            tf2=3; % flag variable
            
        end
        % the signal inputs for the function are 1000 samples of current window
        % and smoothed outputs from previous window(which are availabe from
        % 2nd window onwards, here only present and past signals are used
        % so the system is causal)
        [bpm(i),ty,tf,tf2, prev] = bpmgeneration(sig(:,curSegment),srate,step,ha,f,dl,i,bpm,ty,tf,tf2,prev);
        % the output parameter prev of this function from one window is used
        % as input for the next window
    end
    %**********************************************************************
    
    % Codes to save results
    if(idnb<14)
        plot(BPM0',bpm,'o')
        xlabel('Ground?Truth of Heart Rate (BPM)');
        ylabel('Estimates of Heart Rate (BPM)')
        hold on;
        md(idnb)=mean(abs(BPM0'-bpm));
        mdp(idnb)=mean(abs(BPM0'-bpm)./BPM0')*100;
        
    end
    if (idnb==1)
        y=bpm;
        x=BPM0';
        
    elseif(idnb<14)
        y=[y bpm];
        x=[x BPM0'];
    end
    
    % Codes to clear variables
    clear BPM;
    clear bpm;
    
end
X=x-mean(x);
Y=y-mean(y);
kpc=sum(X.*Y)/sqrt(sum(X.^2)*sum(Y.^2))
md
mdp








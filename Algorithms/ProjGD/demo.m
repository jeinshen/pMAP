% Run install_mex in the PROPACK folder to install first.
% This demo will demonstrate the usages of ProjGD_1D, ProjGD_2D, ProjGD_3D. 

clear; close all; clc;
addpath('./PROPACK/');

prompt = 'Select dimension of the test signal, from 1-3: ';
EX = input(prompt);

% EX=1; % test a 1D example
% EX=2; % test a 2D example
% EX=3; % test a 3D example

switch EX
    
    case 1
        
        % 1D Example
        
        disp('Testing a 1D example...')
        
        % Generate an undamped 1D spectraly sparse signal ox1 of size 1024 with model order 5. There are 
        % enforced separations between the frequencies. It's observed at 128 loacations. Indices for the 
        % samples are stored in K1. The frequencies are stored in f1.
        [K1,ox1,f1] = generate_signal(1024,5,128,'true','false');
        obs1 = ox1(K1); % observed samples stored as a vector
        
        % The maximum number of iteration is set to be 9999. The PGD algorithm will terminate when the 
        % relative change in signal is less than 1e-6 or the gradient magnitude is less than 1e-3. The
        % projection parameter is computed from initialization. The stepsize is chosen via line search.
        % The relative change in signal (relative.change) and the gradient magnitude (gradient.magnitude)
        % per iteration will be displayed.
        [si1,iter1,x1,ratio1,fv1,gm1,step1,t1] = ProjGD_1D(obs1,1024,5,K1,9999,1e-6,1e-3,0,0,1);
        
        % si1 is 0 if failed and 1 if convergence is achieved; iter1 is the number of iteration; x1 is 
        % the recovered signal; ratio1 stores the relative change in signal per iteration; fv1 stores 
        % the function value per iteration; gm1 stores the gradient magnitude per iteration; step1 
        % stores the step size per iteration; t1 stores the time elapsed per iteration.
      
    case 2
        
        % 2D Example
        
        disp('Testing a 2D example...')
        
        % Generate an undamped 2D spectraly sparse signal ox2 of size 64*64 with model order 5. There 
        % are enforced separations between the frequencies. It's observed at 256 loacations. Indices 
        % for the samples are stored in K2. The frequencies for each dimension are stored in the 
        % corresponding column of f2.
        [K2,ox2,f2] = generate_signal([64 64],5,256,'true','false');
        obs2 = ox2(K2); % observed samples stored as a vector
        
        % The maximum number of iteration is set to be 9999. The PGD algorithm will terminate when the 
        % relative change in signal is less than 1e-6 or the gradient magnitude is less than 1e-3. The
        % projection parameter is computed from initialization. The stepsize is chosen by minimizing a 
        % fourth order polynomial per iteration. The relative change in signal (relative.change) and the 
        % gradient magnitude (gradient.magnitude) per iteration will be displayed.
        [si2,iter2,x2,ratio2,fv2,gm2,step2,t2] = ProjGD_2D(obs2,64,64,5,K2,9999,1e-6,1e-3,0,1,1);
        
        % si2 is 0 if failed and 1 if convergence is achieved; iter2 is the number of iteration; x2 is 
        % the recovered signal stored as a vector; ratio2 stores the relative change in signal per 
        % iteration; fv2 stores the function value per iteration; gm2 stores the gradient magnitude per 
        % iteration; step2 stores the step size per iteration; t2 stores the time elapsed per iteration.
        
    case 3
        
        % 3D Example
        
        disp('Testing a 3D example...')
        
        % Generate a damped 3D spectraly sparse signal ox3 of size 64*32*32 with model order 5. There 
        % are no enforced separations between the frequencies. We have 1 percent of all the samples. 
        % Indices for the samples are stored in K3. The frequencies for each dimension are stored in 
        % the corresponding column of f3.
        [K3,ox3,f3] = generate_signal([64 32 32],5,0.01,'false','true');
        obs3 = ox3(K3); % observed samples stored as a vector
        
        % The maximum number of iteration is set to be 9999. The PGD algorithm will terminate when the 
        % relative change in signal is less than 1e-6 or gradient magnitude is less than 1e-3. The 
        % projection parameter is computed from initialization. The stepsize is chosen via line search. 
        % The relative change in signal (relative.change) and the gradient magnitude (gradient.magnitude)
        % per iteration will be displayed.
        [si3,iter3,x3,ratio3,fv3,gm3,step3,t3] = ProjGD_3D(obs3,64,32,32,5,K3,9999,1e-6,1e-3,0,0,1);
        
        % si3 is 0 if failed and 1 if convergence is achieved; iter3 is the number of iteration; x3 is 
        % the recovered signal stored as a vector; ratio3 stores the relative change in signal per 
        % iteration; fv3 stores the function value per iteration; gm3 stores the gradient magnitude per 
        % iteration; step3 stores the step size per iteration; t3 stores the time elapsed per iteration.
        
    otherwise
        
        error('Current supports for dimensions 1-3!')

end
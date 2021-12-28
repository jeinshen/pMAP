% This demo will demonstrate the usages of FIHT_1D, FIHT_2D, FIHT_3D. The
% usages of IHT_1D, IHT_2D, IHT_3D are similar.

clear all;close all;clc;
addpath('./PROPACK/'); % run install_mex to install first

EX=1; % Test a 1D example
%EX=2; % Test a 2D example
%EX=3; % Test a 3D example

switch EX
    
    case 1
        
        % 1D Example
        
        disp('Testing a 1D example...')
        
        % Generate 1D signal ox1 of length 1024 with model order 5 and it's observed
        % at 128 different loacations. The indices for the samples are stored in K1.
        % The frequencies are stored in f1.
        [K1,ox1,f1]=generate_signal(1024,5,128);
        obs1=ox1(K1); % observed samples, stored as a vector
        
        % Run FIHT algorithm, the maximum number of iterations are set to be 500
        % and the algorithm will terminate when the relative error
        % || x_{l}-x{l-1} || / || x_{l-1} || falls below 1e-5.
        tic;[si1,iter1,ratio1,x1]=FIHT_1D(obs1,1024,5,K1,500,1e-5,1);toc;
        
        % si1 is 0 if failed and 1 if successful, iter1 is the number of
        % iterations taken, ratio1 stores the relative mean square error (RMSE)
        % every iteration, x1 is the reconstructed signal.
        
        if si1==1
            plot(1:iter1,log10(ratio1));grid on;
            xlabel('ITER','FontSize',12)
            ylabel('RMSE','FontSize',12)
            title('$(N,r,m)=(1024,5,128)$','Interpreter','latex','FontSize',15)
        else
            disp('Failed!')
        end
        
    case 2
        
        % 2D Example
        
        disp('Testing a 2D example...')
        
        % Generate 2D signal ox2 of dimensions 64*32 with model order 5 and it's
        % observed at 256 different loacations. The vector indices for the samples
        % are stored in K2. The frequencies for each dimension are stored in the
        % corresponding column of f2.
        [K2,ox2,f2]=generate_signal([64 32],5,256);
        obs2=ox2(K2); % observed samples, stored as a vector
        
        % Run FIHT algorithm, the maximum number of iterations are set to be 500
        % and the algorithm will terminate when the relative error
        % || x_{l}-x{l-1} || / || x_{l-1} || falls below 1e-5. The relative error
        % of each iteration will be displayed.
        tic;[si2,iter2,ratio2,x2]=FIHT_2D(obs2,64,32,5,K2,500,1e-5,1);toc;
        
        % si2 is 0 if failed and 1 if successful, iter2 is the number of
        % iterations taken, ratio2 stores the relative mean square error (RMSE)
        % every iteration, x2 is the reconstructed signal.
        
        if si2==1
            plot(1:iter2,log10(ratio2));grid on;
            xlabel('ITER','FontSize',12)
            ylabel('RMSE','FontSize',12)
            title('$(N,r,m)=(64 \times 32,5,256)$','Interpreter','latex','FontSize',15)
        else
            disp('Failed!')
        end
        
    case 3
        
        % 3D Example
        
        disp('Testing a 3D example...')
        
        % Generate 3D signal ox3 of dimensions 16*16*16 of model order 5 and it's
        % observed at 512 different loacations. The vector indices are stored in K3.
        % The frequencies for each dimension are stored in the corresponding column
        % of f3.
        [K3,ox3,f3]=generate_signal([32 32 16],5,1024);
        obs3=ox3(K3); % observed samples, stored as a vector
        
        % Run FIHT algorithm, the maximum number of iterations are set to be 500
        % and the algorithm will terminate when the relative error
        % || x_{l}-x{l-1} || / || x_{l-1} || falls below 1e-5. The relative error
        % of each iteration will be displayed.
        tic;[si3,iter3,ratio3,x3]=FIHT_3D(obs3,32,32,16,5,K3,500,1e-5,1);toc;
        
        % si3 is 0 if failed and 1 if successful, iter3 is the number of
        % iterations taken, ratio3 stores the relative mean square error (RMSE)
        % every iteration, x3 is the reconstructed signal.
        
        if si3==1
            plot(1:iter3,log10(ratio3));grid on;
            xlabel('ITER','FontSize',12)
            ylabel('RMSE','FontSize',12)
            title('$(N,r,m)=(32 \times 32 \times 16,5,1024)$','Interpreter','latex','FontSize',15)
        else
            disp('Failed!')
        end

end
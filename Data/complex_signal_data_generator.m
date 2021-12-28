
function Y = complex_signal_data_generator(N, m, r, theta, instance_amount)


%%
% Generate sample data (spectral spare signals) for Experiment 2
% The data is generated via generate_signal method, which is provided in 
%  the FIHT solver package
%   Case test_data 1 create signal data used for Test 2.a and Test 2.c
%   where no noises is considered
%   Case test_data 2 create signal data for Test 2.b where some
%   observations are with noises
%%

K = zeros(m, instance_amount);
OX = zeros(N, instance_amount);

if theta == 0
    for i = 1 : instance_amount
        [K1,signal,~] = generate_signal(N,r,m,'false','false');
        K(:, i) = K1;
        OX(:, i)= signal;
    end
    
    Y{1} = K; Y{2} = OX;
else
    OX_N = zeros(N, instance_amount);

    for i = 1 : instance_amount
            
        [K1,signal,~] = generate_signal(N,r,m,'false','false');
        mpartial = floor(1 * m /3);
        h=randn(mpartial,1)+1i*randn(mpartial,1);
        e = theta * norm(signal(K1(1:mpartial))) * h / (norm(h));
        e = [e; zeros(m - mpartial, 1)];
        signal_N = signal;
        signal_N(K1) = signal(K1) + e;
        OX_N(:, i) = signal_N;
            
        K(:, i) = K1;
            
        OX(:, i)= signal;
            
    end
        
	Y{1} = K; 
    Y{2} = OX;
	Y{3} = OX_N;

end

% 
% 
% test_data = 2;
% 
% switch test_data
%     
%     case 1
%         
%         r = 120; % Objective rank
%         m = 1200;
%         N = 1999;
%         
%         K = zeros(m, 50);
%         
%         for i = 1 : 50
%             
%             [K1,signal,~] = generate_signal(N,r,m,'false','false');
%             K(:, i) = K1;
%             OX(:, i)= signal;
%             
%         end
%         
%         Y{1} = K; Y{2} = OX;
% 
%         
%     case 2
%         
%         r = 20;
%         theta = 0.2;
%         m = 1200;
%         N = 2000;
% %         noise = 1;
%         
%         K = zeros(m, 50);
%         
%         for i = 1 : 50
%             
%             [K1,signal,~] = generate_signal(N,r,m,'false','false');
%             mpartial = floor(1 * m /3);
%             h=randn(mpartial,1)+1i*randn(mpartial,1);
%             e = theta * norm(signal(K1(1:mpartial))) * h / (norm(h));
%             e = [e; zeros(m - mpartial, 1)];
%             signal_N = signal;
%             signal_N(K1) = signal(K1) + e;
%             OX_N(:, i) = signal_N;
%             
%             K(:, i) = K1;
%             
%             OX(:, i)= signal;
%             
%         end
%         
%         Y{1} = K; Y{2} = OX;
%         Y{3} = OX_N;
% 
% end
end
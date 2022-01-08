% NOMA
close all;
clear all;clc;
warning('off');
rand('twister',mod(floor(now*8640000),2^31-1));

N          = 1;            % array number of users 
K          =2;            % number of users 
SNR_dB     = 10;     % dBW
%%%%% noise
N0=10^((-174-30) / 10); %-174dBm  
B=10^4; %10MHz
 noise_maxpower_original   = N0*B;            % % W
noise_maxpower_original   = 10^((-50-30) / 10);            % % W
% noise_maxpower   = 1;            % % W
trans_maxpower_all =0; % 
% error = 0.05; 
error_t = 0:0.01:0.14; 
% rate_min_dB   = [1,1.2,1.4,1.6,1.8,2,2.2]  ;   %bit/s/hz
rate_min  = [1.2]  ; 
prob=0.05;

%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('H_d_all');
  Power=zeros(length(error_t),1);
   Rate=zeros(K,length(error_t ),10);
 for i_p     = 1 : length(error_t)
 
   error = error_t (i_p);
   num_loop = 100; 
  Power_t=0;
for loop =1 : num_loop
    outerflag=1;      
    H=H_d_all(1:N,1:K,loop)/sqrt(noise_maxpower_original);
    noise_maxpower=1;

    for k=1:K
        H_error(k)=error*norm(H(:,k),'fro');
    end
%%  For different SNR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  loop |  num_J  |  SNR  |  i  |  trans_SNR | relay_SNR \n');

    trans_maxpower=100 ;%trans_maxpower_all(i_p);
  
    F_ini=randn(N,K)*sqrt(trans_maxpower/(N*K));
    F(:,:,1)=full(F_ini);
  
    num_iterative = 15000;
    for n  = 1 : num_iterative
       %%%%%  Optimize F  %%%%%
        [F_1,power_sub,x,innerflag] = Generate_beamforming_F(N,  K, H, H_error,...
                    prob, noise_maxpower, B, rate_min);
        power(n+1)=power_sub;
        F(:,:,n+1)=F_1;
        if innerflag==0
            outerflag=0;
            break;
        end     
       
        fprintf('   %g  |  %g  |  %g  \n',loop,error_t (i_p), n);
        if abs(power(n+1)-power(n))<10^(-4) 
            break;
        end
        
    end
    if outerflag==0
        break;
    end
     Power_t= Power_t+power_sub;


   
    PPower(loop,i_p)=real(power_sub);

    iteration(loop,i_p)=n;

end
 
 Power(i_p)=Power_t/num_loop;
 end

    save('Power','Power');

a=1;
    

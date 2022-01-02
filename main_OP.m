close all;
clear all;clc;
warning('off');
rand('twister',mod(floor(now*8640000),2^31-1));

N          = 1;            
K          =2;            % number of users 
SNR_dB     = 10;     % dBW
%%%%% noise
N0=10^((-174-30) / 10); %-174dBm  
B=10^4; %10MHz
 noise_maxpower_original   = N0*B;            % % W
noise_maxpower_original   = 10^((-50-30) / 10);            % % W
% noise_maxpower   = 1;            % % W
trans_maxpower_all =0; % 
error = 0.1; 
% error_t = 0:0.01:0.14; 

% rate_min_dB   = [1,1.2,1.4,1.6,1.8,2,2.2]  ;   %bit/s/hz
rate_min  =2; 
 prob=0.05;
prob_set=0.02:0.01:0.16; 
%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('H_d_all');
  Power=zeros(length(prob_set),1);
   Rate=zeros(K,length(prob_set ),10);
 for i_p     = 1 : length(prob_set)
 
prob=prob_set(i_p);
   num_loop = 10; 
  Power_t=0;
for loop =1 : num_loop
    outerflag=1;      
    H=H_d_all(1:N,1:K,loop)/sqrt(noise_maxpower_original);
    noise_maxpower=10;
%     H=H_d_all(1:N,1:K,loop);
%     G=G_r_all(1:M,1:N,1:K,loop);
%     noise_maxpower=noise_maxpower_original;
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
                    F(:,:,n),  prob, noise_maxpower, B, rate_min);
        power(n+1)=power_sub;
        F(:,:,n+1)=F_1;
        if innerflag==0
            outerflag=0;
            break;
        end     
       
        fprintf('   %g  |  %g  |  %g  \n',loop,prob_set(i_p), n);
        if abs(power(n+1)-power(n))<10^(-4) 
            break;
        end
        
    end
    if outerflag==0
        break;
    end
     Power_t= Power_t+power_sub;

 F_temp=F_1;
    flag=ones(100,1);
    for i_loop=1:100
      for k=1:K
          H_error_channel(:,k)=sqrt(H_error(k)^2/2)*(unifrnd(-1,1,N) + sqrt(-1)* unifrnd(-1,1,N));
          H_prac(:,k)=H(:,k)+H_error_channel(:,k);
          y(k,1)=norm((H_prac(:,k))*F_temp(:,k),2)^2;
          z_ini(k,1)=norm((H_prac(:,k))*F_temp(:,k:K),2)^2 ...
               +noise_maxpower;
          RATE(k,i_loop)=B*log2(1+y(k,1)/z_ini(k,1));
          if  RATE(k,i_loop) < rate_min
              flag(i_loop)=0;
          end
      end
    end
    Rate(:,loop,i_p)=sum(RATE,2)/i_loop;
    Rate_ratio(loop,i_p)=sum(flag~=0)/length(flag);
   
    PPower(loop,i_p)=real(power_sub);

    iteration(loop,i_p)=n;

end
 
 Power(i_p)=Power_t/num_loop;
 end

    save('Power','Power');
       save('Rate','Rate');
    save('Rate_ratio','Rate_ratio');

a=1;
    

close all;
clear all;clc;
N          = 1;            % array number of users
K          = 4;         % number of users in each group


%%%%% Large scale path loss
PL_0=10^(-30/10); %dB the channel gain at the reference distance

for k=1:K
   x_user_t(k) =3*rand;y_user_t(k)=3*rand;
% x_user_t(1)=30; y_user_t(1)=0;
% x_user_t(2)=35; y_user_t(2)=-5;
% x_user_t(3)=30; y_user_t(3)=2;
% x_user_t(4)=33; y_user_t(4)=1;
% x_user_t(5)=31; y_user_t(5)=2;
% x_user_t(6)=33; y_user_t(6)=-1;
% x_user_t(7)=31; y_user_t(7)=-2;

% d_BI=sqrt((x_irs-x_bs)^2+(y_irs-y_bs)^2); %m distance from the BS to IRS
% pathloss_BI=sqrt(PL_0*(d_BI)^(-2.2));    % Large-scale pass loss from the BS to the IRS

 x_user(k)=30+3*rand;y_user(k)=3*rand;
% x_user(1)=70; y_user(1)=1;
% x_user(2)=65; y_user(2)=-3;
% x_user(3)=70; y_user(3)=3;
% x_user(4)=73; y_user(4)=1;
% x_user(5)=72; y_user(5)=2;
% x_user(6)=70; y_user(6)=2;
% x_user(7)=70; y_user(7)=-1;

end
for k=1:K
    d_BU(k)=sqrt((x_user_t(k)-x_user(k))^2+(y_user_t(k)-y_user(k))^2);%sqrt(d^2+d_v^2);  %m distance from the t_users to the r_users
    pathloss_BU(k)=sqrt(PL_0*(d_BU(k))^(-4));  % Large-scale pass loss from the t_users to the r_users
end

%% Simulation loop %%%%%%%%%%%%%%%%%%%%00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 500; 

for loop = 1 : num_loop
    T1=cputime;


    for k=1:K

        H_d_temp=sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1)); % small scale pass loss from the BS to the user
        H_d_all(:,k,loop)=pathloss_BU(k)*H_d_temp;

    end 
    
end




save('H_d_all','H_d_all');

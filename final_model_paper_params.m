%close all
clear

% Concentrations in mM
% Rate in 1/hr
% Bacte3rial density in OD
rA=0.7; % Maximal Growth rate of 1A01
YAG=1/6; %Yield for GlcNac for 1A01
YAE=1/9.4; % Accetate addition by 1A01 per OD of growth
KAG=10e-3; % Km for 1A01 growth on glcnac

mu_lag_AM=(0.1)/6.2/0.02; % Maximal pyruvate production rate /OD of 1A01 during lag
mu_lag_AE=0.15/0.1*mu_lag_AM; % Maximal acetate production rate /OD of 1A01 during lag
mu_lag_AG=0.1/0.1*mu_lag_AM; % Maximal glcnac consumption rate /OD of 1A01 during lag

mu_str_AM=1.7; % Maximal pyruvate production rate /OD of 1A01 during stress
mu_str_AE=(5.5-2.9)/3.9*mu_str_AM; % Maximal acetate production rate /OD of 1A01 during stress
mu_str_AG=(3.25-0.93)/3.9*mu_str_AM; % Maximal glcnac consumption rate /OD of 1A01 during stress

mu_lag_AM=mu_str_AM; % Maximal pyruvate production rate /OD of 1A01 during lag
mu_lag_AE=mu_str_AE; % Maximal acetate production rate /OD of 1A01 during lag
mu_lag_AG=mu_str_AG;

deltaA=0.5; % Death Rate of 1A01 under acetate stress and GlcNAc runout

rBE=0.35; % Maximal Growth rate of 3B05 on accetate alone outside accetate stress
rBM=0.2; % Maximal Growth rate of 3B05 on pyruvate alone
rB_plus=0.25; % Maximal Growth rate of 3B05 on acetate+pyruvate during stress
KBM=10e-3; % Km for 3B05 growth on pyruvate
KBE=10e-3; % Km for 3B05 growth on accetate
YBE=1/32; % Accetate removal by 3B05 per OD of growth
YBM=3/2*YBE; %1/yield for pyruvate for 3B01 (1st guess is 20)

EA1=3; %Critical accetate for growth arrest of real 1A01 (between 2.5-3)
EA2=4; % (somewhere between 4 and 4.5)
EB1=2.5; % Switching value for growth substrate of 3B05 (somewhere between 3 and 3.5)
EB2=3.5; % Switching value for growth rate of 3B05 (somewhere between 3.5 and 4)
deltaE=1/5; % Slope of tangent function for switching
b=0.75; % Ratio of biomass from acetate

sigma_min=0; % Minimum value of internal variable
sigma_max=1; % Maximum value of internal variable (basically the number of hours of lag)
sigma_c=0.1; % Critical internal variable value below which there is growth
sigma_lag=0.1; % internal variable recovery parameter - taken to be 10 hours
sigma_str=0.5; % internal variable stress paramter

disp=0;

num_cycles=10; % Number of cycles in each run
dt = 1e-4;
cycle_list=24;%:2:24;  % List of cycle durations to iterate over
for time_num=1:length(cycle_list)
    T=cycle_list(time_num);
    num_steps=ceil(T/dt);
    fold_list=40;%^(cycle_list/24);%logspace(1,3,20);
    for fold_num=1:length(fold_list)
        fold=fold_list(fold_num);
        glcnac_list=5;%logspace(0,2,15);
        for glcnac_count=1:length(glcnac_list)
            rhoA=zeros(num_steps,1);
            rhoB=zeros(num_steps,1);
            sigma=zeros(num_steps,1); % Lag time is stored as a function of accetate exposure
            met=zeros(num_steps,1);
            acc=zeros(num_steps,1);
            glc=zeros(num_steps,1);
            
            count=0;
            for cycle=1:num_cycles
                if cycle>1  % Dilution
                    rhoA(1)=rhoA(end)/fold;
                    rhoB(1)=rhoB(end)/fold;
                    acc(1)=acc(end)/fold;
                    met(1)=met(end)/fold;
                    glc(1)=glc(end)/fold;
                    sigma(1)=sigma(end);
                else % Initial condition of experiment
                    x=1/3;
                    rhoA(1)=x*2.0e-2;
                    rhoB(1)=(1-x)*2.0e-2;
                    met(1)=0;
                    acc(1)=0;
                    glc(1)=0;
                    sigma(1)=sigma_min;
                end
                
                % Addition of glcnac
                glc(1)=glc(1)+glcnac_list(glcnac_count);
                t=0;
                
                % Time evolution
                for i=1:num_steps-1
                    if acc(i)>=EA2 % Death
                        sigma(i+1)=min(sigma(i)+(dt*acc(i)*sigma_str),sigma_max);
                        rhoA(i+1)=rhoA(i)-dt*deltaA*rhoA(i);
                        met(i+1)=met(i);
                        acc(i+1)=acc(i);
                        glc(i+1)=glc(i);
                    elseif acc(i)>=EA1  % Acid stress
                        sigma(i+1)=min(sigma(i)+(dt*acc(i)*sigma_str),sigma_max); % Increase lag time for glcnac to catch up to, up to maximum
                        rhoA(i+1)=rhoA(i);
                        met(i+1)=met(i)+dt*mu_str_AM*rhoA(i)*glc(i)/(glc(i)+KAG);
                        acc(i+1)=acc(i)+dt*mu_str_AE*rhoA(i)*glc(i)/(glc(i)+KAG);
                        glc(i+1)=glc(i)-dt*mu_str_AG*rhoA(i)*glc(i)/(glc(i)+KAG);
                    else
                        if sigma(i)<sigma_c % Growth
                            sigma(i+1)=max(sigma(i)-dt*sigma_lag*glc(i)/(glc(i)+KAG),sigma_min); % Increase lag time for glcnac to catch up to, up to maximum
                            rhoA(i+1)=rhoA(i)+dt*rA*rhoA(i)*glc(i)/(glc(i)+KAG);
                            met(i+1)=met(i);
                            acc(i+1)=acc(i)+dt*rA*rhoA(i)/YAE*glc(i)/(glc(i)+KAG);
                            glc(i+1)=glc(i)-dt*rA*rhoA(i)/YAG*glc(i)/(glc(i)+KAG);
                        else % Lag
                            sigma(i+1)=sigma(i)-dt*sigma_lag*glc(i)/(glc(i)+KAG); % Increase lag time for glcnac to catch up to, up to maximum
                            rhoA(i+1)=rhoA(i);
                            met(i+1)=met(i)+dt*mu_lag_AM*rhoA(i)*glc(i)/(glc(i)+KAG);
                            acc(i+1)=acc(i)+dt*mu_lag_AE*rhoA(i)*glc(i)/(glc(i)+KAG);
                            glc(i+1)=glc(i)-dt*mu_lag_AG*rhoA(i)*glc(i)/(glc(i)+KAG);
                        end
                    end
                    
                    if acc(i)<EB1 %Growth on single substrates
                        rhoB(i+1)=rhoB(i)+dt*rhoB(i)*(rBE*acc(i)/(acc(i)+KBE)+rBM*met(i)/(met(i)+KBM)); %% For low acid, Monod growth on accetate
                        met(i+1)=met(i+1)-dt*rhoB(i)*rBM*met(i)/(met(i)+KBM)/YBM;
                        acc(i+1)=acc(i+1)-dt*rhoB(i)*rBE*acc(i)/(acc(i)+KBE)/YBE;
                    else % Combined growth
                        theta_B2=switch_fun(acc(i),rB_plus,0.0,EB2,1/deltaE); %% Accetate growth is a soft tanh function
                        rhoB(i+1)=rhoB(i)+dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM); %% For low acid, Monod growth on accetate
                        %met(i+1)=met(i+1)-(1-b)*dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM)/YBM;
                        %acc(i+1)=acc(i+1)-b*dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM)/YBE;
                        met(i+1)=met(i+1)-(1-b)*dt*rhoB(i)*theta_B2*acc(i)*met(i)/(acc(i)*met(i)+acc(i)*KBM+met(i)*KBE)/YBM;
                        acc(i+1)=acc(i+1)-b*dt*rhoB(i)*theta_B2*acc(i)*met(i)/(acc(i)*met(i)+acc(i)*KBM+met(i)*KBE)/YBE;
                    end
                    
                    if mod(t,1)<dt
                        time_cycle(count+1)=t;
                        rhoA_cycle(count+1)=rhoA(i);
                        rhoB_cycle(count+1)=rhoB(i);
                        acc_cycle(count+1)=acc(i);
                        glc_cycle(count+1)=glc(i);
                        pyr_cycle(count+1)=met(i);
                        count=count+1;
                    end
                    t=t+dt;
                end
                
                %Adding an entry for time=24hrs
                %                 time_cycle(count+1)=t;
                %                 rhoA_cycle(count+1)=rhoA(i);
                %                 rhoB_cycle(count+1)=rhoB(i);
                %                 acc_cycle(count+1)=acc(i);
                %                 glc_cycle(count+1)=glc(i);
                %                 pyr_cycle(count+1)=met(i);
                %                 count=count+1;
            end
        end
    end
end

%% Plot

% Plot final cycle details
figure(1)
%figure('units','normalized','outerposition',[0 0 1 1])
clf;
subplot(2,2,1)
time=(1:num_steps)*dt;
semilogy(time,rhoA,'LineWidth',2,'Color','#0E72BA');
hold on
semilogy(time,rhoB,'LineWidth',2);
xlabel('Time (in hours)')
title('Species Densities (in OD)')
legend('1A01','3B05')

subplot(2,2,2)
cla;
plot(time,glc,'-','LineWidth',2);
hold on
plot(time,met,'-','LineWidth',2);
plot(time,acc,'-','LineWidth',2);
plot(time,sigma,'-','LineWidth',2);
%ylim([1e-6,1e2])
xlabel('Time (in hours)')
title('Nutrients')
legend('glcNac','Pyruvate','Accetate','glcNAc Counter','Lag Time')

%%
figure(3)
clf;
subplot(2,1,2)
cla;
transp=1;
s=plot(time,glc,'--','LineWidth',5,'Color','blue');
s.Color(4)=transp;
hold on
s=plot(time,acc,':','LineWidth',5,'Color','red');
s.Color(4)=transp;
s=plot(time,met,'-','LineWidth',5,'Color','#7e1e9c');
s.Color(4)=transp;
%plot(time,EA1*ones(length(time),1),'-.','LineWidth',5,'Color','red');
% s=patch([0,0,8,8],[0,12,12,0],'blue','EdgeColor','none');
% alpha(s,0.1);
% s=patch([8,8,13.34,13.34],[0,12,12,0],'green','EdgeColor','none');
% alpha(s,0.1);
% % s=patch([14,14,15.2,15.2],[0,1.2,1.2,0],'green','EdgeColor','green');
% % alpha(s,0.1);
% s=patch([13.34,13.34,18.64,18.64],[0,12,12,0],'yellow','EdgeColor','none');
% alpha(s,0.2);
% s=patch([13.34,13.34,16.77,16.77],[0,12,12,0],'red','EdgeColor','none');
% alpha(s,0.2);
% s=patch([16.77,16.77,24,24],[0,12,12,0],'white','EdgeColor','none');
% alpha(s,0.1);
%xlabel('Time (in hours)')
set(gca,'fontsize', 30);
ax=gca;
k=0.015;
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.

xlim([0,24])
ylim([0,5.5])
xticks(0:4:24);
%title('Nutrients')
[~,h]=legend('GlcNAc','Acetate','Acid-induced Metabolites');
xlabel('Time (in hours)')
set(gca,'fontsize', 36);
%set(gca,'YTickLabel',[]);
xlim([0,24])
%ylim([0,1.1])
xticks(0:4:24);
h(4).LineWidth=3;
h(6).LineWidth=3;
h(1).FontSize=30;
h(2).FontSize=30;
h(3).FontSize=30;
legend('boxoff')

%set(gca, 'color', 'none');

subplot(2,1,1)
cla;
yyaxis left
transp=0.75;
s=semilogy(time,rhoA,'-','LineWidth',5,'Color','black');
s.Color(4)=1;
hold on
s=semilogy(time,rhoB,'--','LineWidth',5,'Color','black');
s.Color(4)=1;
ylim([1e-3,1])
set(gca, 'YScale', 'log')
ax=gca;
set(ax,'YColor','black');
k=0.15;
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.
%set(gca, 'color', 'none');
yticks([1e-4,1e-3,1e-2,1e-1,1])


yyaxis right
%plot(time,sigma_c*ones(length(time),1),'-.','LineWidth',5,'Color','#EE1FEF');
%hold on
plot(time,sigma,'-','LineWidth',5,'Color','#EE1FEF');

[~, h] = legend('1A01','3B05','\sigma_A');
% note that even if you plot(x,y,'.') it's a "line" plot
h(4).LineWidth=3;
h(6).LineWidth=3;
h(1).FontSize=30;
h(2).FontSize=30;
h(3).FontSize=30;
legend('boxoff')
ax=gca;
set(ax,'YColor','magenta');

%yyaxis left
% s=patch([0,0,8,8],[0,12,12,0],'blue','EdgeColor','none');
% alpha(s,0.1);
% s=patch([8,8,13.34,13.34],[0,12,12,0],'green','EdgeColor','none');
% alpha(s,0.1);
% % s=patch([14,14,15.2,15.2],[0,1.2,1.2,0],'green','EdgeColor','green');
% % alpha(s,0.1);
% s=patch([13.34,13.34,18.64,18.64],[0,12,12,0],'yellow','EdgeColor','none');
% alpha(s,0.2);
% s=patch([13.34,13.34,16.77,16.77],[0,12,12,0],'red','EdgeColor','none');
% alpha(s,0.2);
% s=patch([16.77,16.77,24,24],[0,12,12,0],'white','EdgeColor','none');
% alpha(s,0.1);
ylim([0,1.2])
xlim([0,24])
set(gca,'XTickLabel',[]);

%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 30);
k=0.015;
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.

%title('Nutrients')

%%
% Plot summary of cycles
figure(222)
clf;
semilogy(rhoA_cycle,'Linewidth',5,'Color','#EF4445')
hold on
semilogy(rhoB_cycle,'Linewidth',5,'Color','#0E72BA')
for t=1:cycle
    xline(t*T,'--','Linewidth',3)
end
legend('1A01','3B05')
xlabel('Time (in hours)')
ylabel('OD')
xlim([0,120])
ylim([0.8e-2,1])
set(gca,'linewidth',2)
set(gca,'fontsize', 24);
%%
subplot(2,1,2)
cla;
plot(acc_cycle,'LineWidth',2)
hold on
plot(glc_cycle,'LineWidth',2)
plot(pyr_cycle,'LineWidth',2)
for t=1:cycle
    xline(t*T,'--','Linewidth',2)
end
legend('Accetate','GlcNAc','Pyruvate')
xlabel('Time (in hours)')
ylabel('Concentration in mM')
xlim([0,240])
set(gca,'fontsize', 24);

writematrix([time_cycle;rhoA_cycle;rhoB_cycle;acc_cycle;glc_cycle;pyr_cycle],'data11_final.csv')

%%
figure(13)
clf;
subplot(2,2,1)
for t=1:length(time)
    %     if glc(t)>1e-3
    %         %env(t)=5-glc(t);
    %         env(t)=glc(t)+met(t)+acc(t)+pyr(t);
    %     else
    %         %env(t)=-met(t);
    %         env(t)=5+0.1*(2-met(t));
    %         env(t)=glc(t)+met(t)+acc(t)+pyr(t);
    %     end
    env(t)=40-(2*acc(t)+8*glc(t)+3*met(t));
end
s=plot((env(3:end-2))/40,(rhoA(4:end-1)-rhoA(3:end-2))./rhoA(3:end-2)/dt,'Linewidth',5,'Color','#EF4445')
alpha(s,.5)
hold on
plot((env(3:end-2))/40,(rhoB(4:end-1)-rhoB(3:end-2))./rhoB(3:end-2)/dt,'Linewidth',5,'Color','#0E72BA')
hold off
alpha(s,.5)
ylabel('Growth Rate (in 1/hr)')
%xlabel('Eco-Coordinate, s (in mM of C)')
xlabel('Cycle Progression')
legend('1A01','3B05')
xlim([0,1])
set(gca,'linewidth',2)
set(gca,'fontsize', 24);

subplot(2,2,2)
plot(time(3:end-2),env(3:end-2),'Linewidth',4,'Color','#006838')
xlabel('Time')
ylabel('Eco-Coordinate, s (in mM of C)')
set(gca,'linewidth',2)
set(gca,'fontsize', 24);
ylim([0,50])

subplot(2,2,3)
plot(time(1:end-1),(rhoA(2:end)-rhoA(1:end-1))./rhoA(1:end-1)/dt,'Linewidth',5,'Color','#EF4445')
hold on
hold on
plot(time(1:end-1),(rhoB(2:end)-rhoB(1:end-1))./rhoB(1:end-1)/dt,'Linewidth',5,'Color','#0E72BA')
ylabel('Growth Rate')
xlabel('Time')
legend('1A01','3B05')
set(gca,'linewidth',2)
set(gca,'fontsize', 24);

subplot(2,2,4)
%%
figure(22)
plot3(glc,acc,met,'Color','#023020','Linewidth',5)
grid on
xlabel('GlcNAC (in mM)')
xlim([0,5])
ylabel('Acetate (in mM)')
ylim([0,3.5])
zlabel('Pyruvate (in mM)')
zlim([0,3])
set(gca,'fontsize', 24);

%%
%T=16;
if 1==0
    num_steps=48/dt;
    rhoA=zeros(num_steps,1);
    rhoB=zeros(num_steps,1);
    sigma=zeros(num_steps,1); % Lag time is stored as a function of accetate exposure
    met=zeros(num_steps,1);
    acc=zeros(num_steps,1);
    glc=zeros(num_steps,1);
    count=0;
    rhoA(1)=1e-2;
    rhoB(1)=1e-2;
    met(1)=0;
    acc(1)=0;
    glc(1)=glcnac_list;
    sigma(1)=sigma_min;
    
    % Time evolution
    for i=1:num_steps-1
        if acc(i)>=EA2 % Death
            sigma(i+1)=min(sigma(i)+(dt*acc(i)*sigma_str),sigma_max);
            rhoA(i+1)=rhoA(i)-dt*deltaA*rhoA(i);
            met(i+1)=met(i);
            acc(i+1)=acc(i);
            glc(i+1)=glc(i);
        elseif acc(i)>=EA1  % Acid stress
            sigma(i+1)=min(sigma(i)+(dt*acc(i)*sigma_str),sigma_max); % Increase lag time for glcnac to catch up to, up to maximum
            rhoA(i+1)=rhoA(i);
            met(i+1)=met(i)+dt*mu_str_AM*rhoA(i)*glc(i)/(glc(i)+KAG);
            acc(i+1)=acc(i)+dt*mu_str_AE*rhoA(i)*glc(i)/(glc(i)+KAG);
            glc(i+1)=glc(i)-dt*mu_str_AG*rhoA(i)*glc(i)/(glc(i)+KAG);
        else
            if sigma(i)<sigma_c % Growth
                sigma(i+1)=max(sigma(i)-dt*sigma_lag*glc(i)/(glc(i)+KAG),sigma_min); % Increase lag time for glcnac to catch up to, up to maximum
                rhoA(i+1)=rhoA(i)+dt*rA*rhoA(i)*glc(i)/(glc(i)+KAG);
                met(i+1)=met(i);
                acc(i+1)=acc(i)+dt*rA*rhoA(i)/YAE*glc(i)/(glc(i)+KAG);
                glc(i+1)=glc(i)-dt*rA*rhoA(i)/YAG*glc(i)/(glc(i)+KAG);
            else % Lag
                sigma(i+1)=sigma(i)-dt*sigma_lag*glc(i)/(glc(i)+KAG); % Increase lag time for glcnac to catch up to, up to maximum
                rhoA(i+1)=rhoA(i);
                met(i+1)=met(i)+dt*mu_lag_AM*rhoA(i)*glc(i)/(glc(i)+KAG);
                acc(i+1)=acc(i)+dt*mu_lag_AE*rhoA(i)*glc(i)/(glc(i)+KAG);
                glc(i+1)=glc(i)-dt*mu_lag_AG*rhoA(i)*glc(i)/(glc(i)+KAG);
            end
        end
        
        if acc(i)<EB1 %Growth on single substrates
            rhoB(i+1)=rhoB(i)+dt*rhoB(i)*(rBE*acc(i)/(acc(i)+KBE)+rBM*met(i)/(met(i)+KBM)); %% For low acid, Monod growth on accetate
            met(i+1)=met(i+1)-dt*rhoB(i)*rBM*met(i)/(met(i)+KBM)/YBM;
            acc(i+1)=acc(i+1)-dt*rhoB(i)*rBE*acc(i)/(acc(i)+KBE)/YBE;
        else % Combined growth
            theta_B2=switch_fun(acc(i),rB_plus,0.0,EB2,1/deltaE); %% Accetate growth is a soft tanh function
            rhoB(i+1)=rhoB(i)+dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM); %% For low acid, Monod growth on accetate
            %met(i+1)=met(i+1)-(1-b)*dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM)/YBM;
            %acc(i+1)=acc(i+1)-b*dt*rhoB(i)*theta_B2*acc(i)/(acc(i)+KBE)*met(i)/(met(i)+KBM)/YBE;
            met(i+1)=met(i+1)-(1-b)*dt*rhoB(i)*theta_B2*acc(i)*met(i)/(acc(i)*met(i)+acc(i)*KBM+met(i)*KBE)/YBM;
            acc(i+1)=acc(i+1)-b*dt*rhoB(i)*theta_B2*acc(i)*met(i)/(acc(i)*met(i)+acc(i)*KBM+met(i)*KBE)/YBE;
        end
        
        rhoA(i+1)=rhoA(i+1)-dt*log(fold)/T*rhoA(i);
        rhoB(i+1)=rhoB(i+1)-dt*log(fold)/T*rhoB(i);
        met(i+1)=met(i+1)-dt*log(fold)/T*met(i);
        acc(i+1)=acc(i+1)-dt*log(fold)/T*acc(i);
        glc(i+1)=glc(i+1)-dt*log(fold)/T*glc(i)+glcnac_list*dt/T;
        
    end
end
%%
figure(5)
clf;
semilogy((1:num_steps)*dt,rhoA,'LineWidth',4)
hold on
semilogy((1:num_steps)*dt,rhoB,'LineWidth',4)

semilogy((1:num_steps)*dt,met,'LineWidth',4)

semilogy((1:num_steps)*dt,glc,'LineWidth',4)

semilogy((1:num_steps)*dt,acc,'LineWidth',4)

semilogy((1:num_steps)*dt,sigma+1e-3,'LineWidth',4)
legend('\rho_A (in OD)','\rho_B (in OD)','Pyruvate (in mM)','GlcNAc (in mM)','Acetate (in mM)', '\sigma')
ylim([1e-3,glcnac_list])
set(gca,'fontsize', 24);
xlabel('Time (in hours)')

%%
figure(15)
clf;
lw=3;
subplot(3,2,1)
title('Stable Cycle')
semilogy(rhoA_cycle(97:120)*8e8,'k-','LineWidth',lw)
hold on
semilogy(rhoB_cycle(97:120)*1.3e9,'k--','LineWidth',lw)
ylabel('Cell Density (/mL)')
legend('1A01','3B05')
ylim([1e6,1e9])
xlim([0,24])
xticks([0,12,24])
yticks(logspace(6,9,4))
set(gca,'fontsize', 15);

subplot(3,2,3)
plot(time(1:end-1),(rhoA(2:end)-rhoA(1:end-1))./rhoA(1:end-1)/dt,'k-','LineWidth',lw)
hold on
plot(time(1:end-1),(rhoB(2:end)-rhoB(1:end-1))./rhoB(1:end-1)/dt,'k--','LineWidth',lw)
ylabel('Growth Rate (/h)')
%ylim([1e6,1e9])
xlim([0,24])
xticks([0,12,24])
set(gca,'fontsize', 15);



subplot(3,2,2)
title('End of Cycle')
semilogy(rhoA_cycle(24:24:120)*8e8,'k-','LineWidth',lw)
hold on
semilogy(rhoB_cycle(24:24:120)*1.3e9,'k--','LineWidth',lw)
ylabel('Cell Density (/mL)')
ylim([1e6,1e9])
xlim([0,5])
xticks(0:5)
yticks(logspace(6,9,4))
set(gca,'fontsize', 15);

subplot(3,2,4)
plot(log(rhoA_cycle(24:24:120)./rhoA_cycle(1:24:97))/24,'k-','LineWidth',lw)
hold on
plot(log(rhoB_cycle(24:24:120)./rhoB_cycle(1:24:97))/24,'k--','LineWidth',lw)
ylabel('Average Growth Rate (/h)')
%ylim([1e6,1e9])
xticks(0:5)
xlim([0,5])
set(gca,'fontsize', 15);

subplot(3,2,6)
ms=10;
title('Stable Cycle')
plot(glc_cycle(24:24:120),'b--^','LineWidth',lw,'MarkerSize',ms)
hold on
plot(acc_cycle(24:24:120),'r:s','LineWidth',lw,'MarkerSize',ms)
plot(pyr_cycle(24:24:120),'s-','Color','#800080','LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','#800080')
ylabel('Concentrations (mM)')
legend('GlcNAc','Acetate','Acid Ind Met')
ylim([0,5.2])
xlim([0,5])
yticks(0:5)
xticks(0:5)
set(gca,'fontsize', 15);
xlabel('Time (h)')

subplot(3,2,5)
title('Stable Cycle')
plot(glc_cycle(97:120),'b--','LineWidth',lw)
hold on
plot(acc_cycle(97:120),'r:','LineWidth',lw)
plot(pyr_cycle(97:120),'Color','#800080','LineWidth',lw)
ylabel('Concentrations (mM)')
legend('GlcNAc','Acetate','Acid Ind Met')
ylim([0,5.2])
xlim([0,24])
yticks(0:5)
xticks([0,12,24])
set(gca,'fontsize', 15);
xlabel('Time (h)')

%%
figure(16)
clf;
lw=3;
subplot(3,3,1)
semilogy(0:23,rhoA_cycle(1:24)*8e8,'k-','LineWidth',lw)
hold on
semilogy(0:23,rhoB_cycle(1:24)*1.3e9,'k--','LineWidth',lw)
ylabel('Cell Density (/mL)')
ylim([1e4,1e9])
xlim([0,24])
xticks([0,12,24])
yticks(logspace(4,9,6))
set(gca,'fontsize', 15);

subplot(3,3,2)
semilogy(0:23,rhoA_cycle(25:48)*8e8,'k-','LineWidth',lw)
hold on
semilogy(0:23,rhoB_cycle(25:48)*1.3e9,'k--','LineWidth',lw)
ylabel('Cell Density (/mL)')
ylim([1e4,1e9])
xlim([0,24])
xticks([0,12,24])
yticks(logspace(4,9,6))
set(gca,'fontsize', 15);

subplot(3,3,3)
semilogy(0:23,rhoA_cycle(25:48)*8e8,'k-','LineWidth',lw)
hold on
semilogy(0:23,rhoB_cycle(25:48)*1.3e9,'k--','LineWidth',lw)
ylabel('Cell Density (/mL)')
ylim([1e4,1e9])
xlim([0,24])
xticks([0,12,24])
yticks(logspace(4,9,6))
set(gca,'fontsize', 15);

subplot(3,3,4)
plot(0:22,log(rhoA_cycle(2:24)./rhoA_cycle(1:23)),'k-','LineWidth',lw)
hold on
plot(0:22,log(rhoB_cycle(2:24)./rhoB_cycle(1:23)),'k--','LineWidth',lw)
ylabel('Growth/Death Rate (/h)')
%ylim([1e6,1e9])
xlim([0,24])
xticks([0,12,24])
set(gca,'fontsize', 15);

subplot(3,3,5)
plot(0:22,log(rhoA_cycle(2+24:24+24)./rhoA_cycle(1+24:23+24)),'k-','LineWidth',lw)
hold on
plot(0:22,log(rhoB_cycle(2+24:24+24)./rhoB_cycle(1+24:23+24)),'k--','LineWidth',lw)
ylabel('Growth Rate (/h)')
%ylim([1e6,1e9])
xlim([0,24])
xticks([0,12,24])
set(gca,'fontsize', 15);

subplot(3,3,6)
plot(0:22,log(rhoA_cycle(2+48:24+48)./rhoA_cycle(1+48:23+48)),'k-','LineWidth',lw)
hold on
plot(0:22,log(rhoB_cycle(2+48:24+48)./rhoB_cycle(1+48:23+48)),'k--','LineWidth',lw)
ylabel('Growth Rate (/h)')
%ylim([1e6,1e9])
xlim([0,24])
xticks([0,12,24])
set(gca,'fontsize', 15);

subplot(3,3,7)
title('Stable Cycle')
plot(0:23,glc_cycle(1:24),'b--','LineWidth',lw)
hold on
plot(0:23,acc_cycle(1:24),'r:','LineWidth',lw)
plot(0:23,pyr_cycle(1:24),'Color','#800080','LineWidth',lw)
ylabel('Concentrations (mM)')
ylim([0,5.2])
xlim([0,24])
yticks(0:5)
xticks([0,12,24])
set(gca,'fontsize', 15);
xlabel('Time (h)')

subplot(3,3,8)
title('Stable Cycle')
plot(0:23,glc_cycle(25:48),'b--','LineWidth',lw)
hold on
plot(0:23,acc_cycle(25:48),'r:','LineWidth',lw)
plot(0:23,pyr_cycle(25:48),'Color','#800080','LineWidth',lw)
ylabel('Concentrations (mM)')
ylim([0,5.2])
xlim([0,24])
yticks(0:5)
xticks([0,12,24])
set(gca,'fontsize', 15);
xlabel('Time (h)')

subplot(3,3,9)
title('Stable Cycle')
plot(0:23,glc_cycle(49:72),'b--','LineWidth',lw)
hold on
plot(0:23,acc_cycle(49:72),'r:','LineWidth',lw)
plot(0:23,pyr_cycle(49:72),'Color','#800080','LineWidth',lw)
ylabel('Concentrations (mM)')
ylim([0,5.2])
xlim([0,24])
yticks(0:5)
xticks([0,12,24])
set(gca,'fontsize', 15);
xlabel('Time (h)')

%%
figure(39)
clf;
plot(time(1:end-1),(rhoA(2:end)-rhoA(1:end-1))./rhoA(1:end-1)/dt,'k-','Linewidth',5)
hold on
plot(time(1:end-1),(rhoB(2:end)-rhoB(1:end-1))./rhoB(1:end-1)/dt,'k--','Linewidth',5)
ylabel('Growth Rate (1/h)')
xlabel('Time (h)')
%legend('1A01','3B05')
xticks([0,6,12,18,24])
set(gca,'linewidth',2)
set(gcf, 'Color', 'None')
set(gca,'fontsize', 24);
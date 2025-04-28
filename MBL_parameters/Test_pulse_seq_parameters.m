clear
close all

J = 1;%32.76*10^3; % krad/s J(F-F)
h = 1;%6.13*10^3; %krad/s 

tau = 7.5; %us, common tau value

% g_vals = [0,0.02,0.16,0.02,0.06,0.14];%.*(1/h);
% v_vals = [0,0,0.004,0.01,0.01,0.01];%.*(1/J);

a_vals = linspace(0,0.2,100); % x disorder strength
b_vals= 0; % y disorder strength
c_vals = linspace(0,0.2,100); % z disorder strength

% a_vals = linspace(0,0.2,100); % x disorder strength
% b_vals= 0; % y disorder strength
% c_vals = linspace(0,0.2,100); % z disorder strength

%v_vals = linspace(0,sqrt(0.25/J^2),100);%[0:0.1:0.5]; % ZZ interaction strength
v_vals = linspace(0,0.2,100);%[0:0.1:0.5]; % ZZ interaction strength
w_vals = -v_vals;
% u_vals = sqrt((0.25/J)-v_vals.^2); % exchange interaction coefficients (u+v)/2, (v-u)/2

% T1 = T*(1+3C-v+w);
% T2 = T*(1+3b-u+v);
% T3 = T*(1- 3a+u-w),
% T11 =T*(1  - 3c - v + w),
% T22 = T*( - 3b - u + v)
% T33 = T*(1+ 3a+u-w),

%% Need to test if any of the Ti are negative  given the above values of the parameters 
flag_neg = 0; 
flag_nonneg = 0;

%% Need to test if any of the Ti are negative  given the above values of the parameters 
for v_index = 1:length(v_vals)
    v = v_vals(v_index);
    w = -v;
    u=0.25;
    %u = sqrt(0.25-((J*v)^2))*(1/J);
        
    for c_index = 1:length(c_vals)
            c = c_vals(c_index);
            b=0;
            a=c;
                 %calculate time intervals 
             T(v_index, c_index, 1) = tau*(1 + 3*c - v + w);
             T(v_index, c_index, 2) = tau*(1+ 3*b - u + v);
             T(v_index, c_index, 3) = tau*(1- 3*a + u - w);
             T(v_index, c_index, 4) = tau*(1 - 3*c - v + w);
             T(v_index, c_index, 5) = tau*(1 - 3*b - u + v);
             T(v_index, c_index, 6) = tau*(1 + 3*a + u - w);

    end 
end

% flag_neg
% flag_nonneg

v_vals_scaled = v_vals .* J; 
c_vals_scaled = c_vals .* h;

% for i=1:6
%     subplot(2,3,i)
%     imagesc(v_vals,c_vals,T(:,:,i))
%     colorbar
%     xlabel('v vals')
%     ylabel('g vals')
% end
% T_positive = T; 
% T_positive(T_positive<=0)=0;

% flag
% for i=1:6
%     figure(1)
%     subplot(2,3,i)
%     imagesc(v_vals,c_vals,(reshape(T(:,2,2,:,i),[length(c_vals),length(v_vals)]))');
%     colorbar
%     set(gca,'YDir','normal')
%     colormap('hot');
%     %clim([0,5])
%     xlabel('v vals')
%     ylabel('g vals')
% 
%     figure(2)
%     subplot(2,3,i)
%     imagesc(v_vals_scaled,c_vals_scaled,(reshape(T(:,2,2,:,i),[length(c_vals),length(v_vals)]))');
%     colorbar
%     set(gca,'YDir','normal')
%     colormap('hot');
%     xlabel('v vals*J')
%     ylabel('g vals*h')
% end

for i=1:6
    figure(1)
    subplot(2,3,i)
    imagesc(v_vals,c_vals,(reshape(T(:,:,i),[length(c_vals),length(v_vals)]))');
    % pcolor(v_vals,c_vals,(reshape(T(:,:,i),[length(c_vals),length(v_vals)]))');
    colorbar
    set(gca,'YDir','normal')
    colormap('hot');
    %clim([0,5])
    xlabel('v')
    ylabel('g')

    figure(2)
    subplot(2,3,i)
    imagesc(v_vals_scaled,c_vals_scaled,(reshape(T(:,:,i),[length(c_vals),length(v_vals)]))');
    % pcolor(v_vals_scaled,c_vals_scaled,(reshape(T(:,:,i),[length(c_vals),length(v_vals)]))');
    colorbar
    set(gca,'YDir','normal')
    colormap('hot');
    xlabel('v vals*J')
    ylabel('g vals*h')
end

titlenames = {'\tau_1','\tau_2','\tau_3','\tau_4','\tau_5','\tau_6'};
figure(1)
for i=1:6
    subplot(2,3,i)
    title(titlenames{i});
    set(findall(gca,'-property','FontSize'),'FontSize',15,'FontName','Times');
    box on
    set(gca,'LineWidth',1)
    %clim([3,15])
    ax = gca;
    ax.XAxis.Exponent = -2;
end


% cycle_time = sum(T,3);
% figure
% imagesc(v_vals,c_vals,cycle_time)
% colorbar
% xlabel('v vals')
% ylabel('g vals')

% figure
% imagesc(v_vals,c_vals,(reshape(sign_indicator(:,1,1,:),[length(c_vals),length(v_vals)]))');
% colorbar
% set(gca,'YDir','normal')
% xlabel('v vals')
% ylabel('g vals')
% 
% figure
% imagesc(v_vals_scaled,c_vals_scaled,(reshape(sign_indicator(:,1,1,:),[length(c_vals),length(v_vals)]))');
% colorbar
% set(gca,'YDir','normal')
% xlabel('v vals*J')
% ylabel('g vals*h')
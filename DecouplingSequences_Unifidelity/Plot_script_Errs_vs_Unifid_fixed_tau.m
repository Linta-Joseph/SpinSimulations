clear
close all

p1_index = 1; 
coupling_index = [1,2,3,4];
delta_index = 1;
overrot_index = 1;
tau_index = 1; 
phasetrans_index = 1; 

filename =['2024-11-08_GaussianSampling(0,1)_8spins_16couplings__P1_0_Tau_4us_' ...
    '4_Coupmaxvals_100_sample_logspacing_AVG_DELTA_onlyBR24.mat'];

load(filename)

styles = {'-o','-s', '-d','->','-p','-*','-x','-h'};
markersize=[7,9,11,13,15];
newcolors = [ 0       0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];


%% image plots of Unitary_fidelity_faircomp_couplingavg for pairs of errors
% order -- sequence, p1, tau, delta, phasetrans,overrot

for i=1:length(sequence_list)      
    colororder(newcolors)

for j=1:length(Couplingmax_vals)

    coupling_index_val=coupling_index(j);
    %subplot(2,4,i)  

    y_var1=reshape(abs(Unitary_fidelity_faircomp_couplingandsampleavg ...
        (i,p1_index,tau_index,coupling_index_val,:)),[1,length(Delta_vals)]);%reshaping the y variable
   
    %plot(Overrot_vals,1-y_var1,styles{i},'LineWidth',2,'MarkerSize',7,'DisplayName',['\tau=',num2str(Tau_vals(tau_index_val)*10^6),'\mus']);%['Coupmax=',num2str(Couplingmax_vals(coupling_index_val)/(2*pi*10^3)),'KHz']);
    plot(Delta_vals/(2*pi),1-y_var1,styles{i},'LineWidth',1,'MarkerSize',7,'DisplayName', sequence_list{i});
    %'DisplayName',['C=',num2str(Couplingmax_vals(coupling_index_val)/(2*pi*10^3)),'KHz']); % [,'MarkerSize',7);%%,['\epsilon=',num2str(Overrot_vals(overrot_index_val))]);%*10^6),'\mus']);%);
    
    hold on 
    clear y_var1

    %if(i==1)
        xlabel('\Delta (Hz)')
        ylabel('1 - F')
    %end
end

end

% Plot formatting
for i=5:5 %variables
    figure(i)
    set(gcf,'Units','centimeters');
    set(gcf,'paperposition',[0 0 4 4],'papersize', [4 4])

    for j=1:1
        %subplot(2,4,j)
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

        box on
        set(gca,'LineWidth',1)
        set(findall(gca,'-property','FontSize'),'FontSize',12);%,'FontWeight','Bold');
        %if (j==4)
             legend('NumColumns',2,'FontSize',12,'Location','best','EdgeColor','None')
        %end

    end
end

   


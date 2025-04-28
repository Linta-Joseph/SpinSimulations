clear
close all 

%% PLOTS
styles = {'-o','-s', '-d','->','-p','-*','-x','-h','-.','-o'};

load('Magnusterms_Iterative.mat')

for n =1:length(sequence_list)
    sequenceName = sequence_list{n};
    testVarLabel = '\tau(\mu s)';

    figure(1)
    %subplot(3,2,1)
    hold on
    plot(testVars,results_norm_hsizes_DIP(n,:,1),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip}^{0}')
    
    figure(2)
    hold on
    %subplot(3,2,2)
    plot(testVars,results_norm_hsizes_DIP(n,:,2),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip}^{1}')
    
    figure(3)
    hold on
    %subplot(3,2,2)
    plot(testVars,results_norm_hsizes_DIP(n,:,3),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip}^{2}')
    
    figure(4)
    hold on
    %subplot(3,2,2)
    plot(testVars,results_norm_hsizes_DIP(n,:,4),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip}^{3}')
    
    figure(5)
    hold on
    %subplot(3,2,2)
    plot(testVars,results_norm_hsizes_DIP(n,:,5),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip}^{4}')

    % figure(6)
    % hold on
    % %subplot(3,2,2)
    % plot(testVars,results_norm_hsizes_DIP(n,:,6),styles{n},...
    %     'DisplayName',sequenceName,LineWidth=3)
    % xlabel(testVarLabel)
    % ylabel('H_{dip}^{5}')
    
    figure(7)
    %subplot(3,2,3)
    plot(testVars,results_norm_hsizes_DELTA(n,:,1),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    hold on
    xlabel(testVarLabel)
    ylabel('H_{\Delta}^{0}')
    
    figure(8)
    %subplot(3,2,5)
    hold on
    plot(testVars,results_norm_hsizes_DELTA(n,:,2),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{\Delta}^{1}')
    
    figure(9)
    %subplot(3,2,4)
    plot(testVars,results_norm_hsizes_CROSS(n,:,2),styles{n},...
        'DisplayName',sequenceName,LineWidth=3)
    xlabel(testVarLabel)
    ylabel('H_{dip-\Delta}^{1}')
    hold on
end 

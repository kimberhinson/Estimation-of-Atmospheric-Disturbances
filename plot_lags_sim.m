clear all
close all

for j = 1:50
    k=0;
i_vect = [150];
for i = i_vect
k = k+1;
clearvars Q R

load(['.\Results\SIM\M_lags' num2str(i) '_' num2str(j)],'R','Q')

for ii = 1:4
    for hh = 1:4
    eval(['Q_vec' num2str(ii) num2str(hh) '{j}(k) = Q(' num2str(ii) ',' num2str(hh) ');']);
    end
end
for ii = 1:3
    for hh = 1:3
    eval(['R_vec' num2str(ii) num2str(hh) '{j}(k) = R(' num2str(ii) ',' num2str(hh) ');']);
    end
end


for ii = 1:4
    for hh = 1:4
        eval(['Q' num2str(ii) num2str(hh) '_end{k}(j) = [Q_vec' num2str(ii) num2str(hh) '{j}(end)];']);         
    end
end

for ii = 1:3
    for hh = 1:3
        eval(['R' num2str(ii) num2str(hh) '_end{k}(j) = [R_vec' num2str(ii) num2str(hh) '{j}(end)];']); 
    end
end
end


end 

meanRb = [];
stdR = [];
for i = 1:3
    meanRb(i) = mean(eval(['R' num2str(i) num2str(i) '_end{1}']));
    stdR(i) = std(eval(['R'  num2str(i) num2str(i) '_end{1}']));   
end

trueQ = 0.01*diag([1,2,6,3]);                                                              
trueR = 0.00001*diag([9,6,4]);

ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};

figure(1)
hold on; grid on
plot([1:3],meanRb,"-s","MarkerSize",5,"MarkerEdgeColor","[0 0.4470 0.7410]","MarkerFaceColor",[0.65 0.85 0.90],"LineStyle","none");
errorbar([1:3],meanRb,-stdR,stdR,"LineStyle","none","Color","[0 0.4470 0.7410]")
hold on; grid on;
xlabel('Diagonal Elements of R')
plot([1:3],diag(trueR),'xk');
legend('Mean of solutions','Standard deviation of solutions','Truth')  
xlim([0 4])
xticks([1:1:3])
ylim([0 15e-5])

meanQb = [];
stdQ = [];
for i = 1:4
    meanQb(i) = mean(eval(['Q' num2str(i) num2str(i) '_end{1}']));
    stdQ(i) = std(eval(['Q'  num2str(i) num2str(i) '_end{1}']));   
end


ccl = {[240/255,248/255,255/255],[240/255,255/255,240/255],[255/255,240/255,245/255],[230/255,230/255,250/255]};
ccd = {[0 0.4470 0.7410],[34/255,139/255,34/255],[139/255 0 0],[75/255,0,130/255]};
    
figure(2)
hold  on; grid on
plot([1:4],meanQb,"-s","MarkerSize",5,"MarkerEdgeColor","[0 0.4470 0.7410]","MarkerFaceColor",[0.65 0.85 0.90],"LineStyle","none");
errorbar([1:4],meanQb,-stdQ,stdQ,"LineStyle","none","Color","[0 0.4470 0.7410]")
hold on; grid on;
xlabel('Diagonal Elements of Q')
plot([1:4],diag(trueQ),'xk');
legend('Mean of solutions','Standard deviation of solutions','Truth') 
xlim([0 5])
xticks([1:1:4])

save(['./Results/SIM/meanQR'],'meanQb','meanRb','stdQ','stdR');

%% .png creation

% figure(1)
% set(gcf,'Units','inches');
% set(gcf,'Position',[0 1 12 5 ])
% screenposition = get(gcf,'Position');
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpng -painters -r400 R_sim_legend
% 
% figure(2)
% set(gcf,'Units','inches');
% set(gcf,'Position',[0 1 12 5 ])
% screenposition = get(gcf,'Position');
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpng -painters -r400 Q_sim_legend

clear all;

k=0;
i_vect = [150];
for i = i_vect

k = k+1;
clearvars Q R A U V S
load(['.\Results\SIM\M_lags' num2str(i) '_1.mat']);

[U,S,V] = svd(A);
S_diag(:,k) = diag(S);
cc(k) = cond(A);
end


figure
histogram(log10(S_diag),'BinWidth',0.5)
hold on; grid on
xlabel('log$\sigma_{\mathcal{A}}$','Interpreter','LaTeX')
ylim([0 5])
%% .png creation

% set(gcf,'Units','inches');
% set(gcf,'Position',[0 1 12 5 ])
% screenposition = get(gcf,'Position');
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpng -painters -r400 Asv_sim
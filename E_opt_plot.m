NG = norm_diff_ALS;
NGp = norm_diff_ALSp;
NGg = norm_diff_ALSg;

figure
for i = 1

semilogy(i,NGp(i,1),'Marker','o','LineStyle','none','Markersize',10,'Color',[0 0.4470 0.7410])
hold on; grid on;
semilogy(i,NG(i,1),'Marker','s','LineStyle','none','Markersize',10,'Color',[0.8500 0.3250 0.0980])
semilogy(i,NGg(i,1),'Marker','p','Color',[0.4660, 0.6740, 0.1880],'LineStyle','none','Markersize',10,'Color',[0.4660, 0.6740, 0.1880])


end

ylabel('E_{opt}')

xtickangle(40)
legend('Ideal', 'ALS', 'ALS colored')

% set(gcf,'Units','inches');
% set(gcf,'Position',[0 1 6 4])
% screenposition = get(gcf,'Position');
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpng -painters -r400 Eopt_sim_compare

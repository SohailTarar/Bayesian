function plot = plot_learning_rate(parameter)
parameter_mean = mean(parameter);
[M,edges] = histcounts(parameter);
[~, id] = max(M);
parameter_mode = edges(id);
d = sort(parameter);
parameter_median = d(end/2);
phi_CI = HDI(parameter,0.95);
histogram(parameter,5000,'normalization','pdf','EdgeColor','none'); xlim([0.93 1.02]); 
hold on; xl = xlim(); yl = ylim(); xWidth = xl(2)-xl(1); yHeight = yl(2)-yl(1);
x = (xl(1) + 0.8 * xWidth); y = (yl(1) + 0.8 * yHeight);
data = sprintf('Mode = %g \n%s %g',parameter_mode,'Median = ',parameter_median);
text(x,y,data); 
lower = phi_CI(1,1); upper = phi_CI(1,2);
yL = get(gca,'YLim');
line([parameter_mean parameter_mean],yL,'LineWidth',2,'Color','red','LineStyle','--'); 
line([lower lower],yL,'LineWidth',2,'Color','red','LineStyle','--'); 
line([upper upper],yL,'LineWidth',2,'Color','red','LineStyle','--'); 
text(parameter_mean,y*1.2,sprintf('Mean = %g',parameter_mean),'HorizontalAlignment','center')
text(lower,y*1.1,sprintf('HDI low = %g',lower),'HorizontalAlignment','right')
text(upper,y*1.1,sprintf('HDI up = %g',upper),'HorizontalAlignment','left')
set(findall(gcf,'-property','FontSize'),'FontSize',16); hold off;
end
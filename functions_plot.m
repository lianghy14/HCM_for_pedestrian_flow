function pfuns = functions_plot
pfuns.potential = @plot_potential;
end


function number = plot_potential(potential,bj_calODH,x,y,zlim,fig_title)
[X,Y] = meshgrid(x,y);
potential(bj_calODH==0) = nan;
plotrange=[min(x),max(x),min(y),max(y),min(min(potential)),zlim];

figure;
subplot(2,1,1);
set(surf(X,Y,potential),'Facecolor','none');
axis(plotrange);
title(fig_title)
xlabel('x(m)'); ylabel('y(m)');
subplot(2,1,2); 
imagesc(x,y, potential,'alphadata',~isnan(potential)); colorbar;
caxis([0 zlim]);
set(gca,'YDir','normal','color',0*[1 1 1]);
xlabel('x(m)'); ylabel('y(m)'); 

number = 0;
end
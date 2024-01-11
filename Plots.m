close all
clear all
clc

blue = ([0 114 178]/255);
bluea = ([0 114 178 180]/255);
reda = ([150 10 10 180]/255);
red = ([150 10 10]/255);

tab=readtable('depth_qual.csv');
models={"FLake","GLM", "GOTM", "Simstrat"};
metrics= {'bias','mae','mse','rmse','ame','mrme','nme', ...
    'r','nse'};
uom= {'ME [°C]','MAE [°C]','MSE [°C^2]','RMSE [°C]','AME [°C]','MRME [-]','NME [-]', ...
    'r [-]','NSE [-]'};
nm=length(models);
nmet=length(metrics);

figure('Units','centimeters','PaperSize',[15,15],'PaperPosition',[0 0 15 15])
ha = tight_subplot(nm,nmet,[.01 .02],[.12 .01],[.095 .01]);

% for j=1:nmet
%     idx=find(strcmp(tab.Metric,metrics{j}));
%     xxplot(j,:)=[min(tab.Value(idx)) max(tab.Value(idx))];
% end
% xxplot=[-2.5 0.5;0.2 2.5; 0 12;0.4 3.5;...
%     0.5 8; -0.2 0.04; -0.25 0.05; 0.85 1; 0 1];
xxplot=[-0.6 0.5;0.2 1; 0 1.5;0.4 1.2;...
    1 5.2; -0.1 0.04; -0.1 0.05; 0.95 1; 0.8 1];
c=0;
for i=1:nm
    for j=1:nmet
        c=c+1;
        idx=find(strcmp(tab.Model,models{i}) & strcmp(tab.Metric,metrics{j}));
        z=tab.Depth(idx);
        val=tab.Value(idx);
        axes(ha(c));
        plot(val,z,'-','color',blue,'linewidth',1.1); set(gca,'ydir','reverse'); hold on
        ii=find(val~=-Inf); val=val(ii);
        if i==2
           disp([metrics{j} ' ' num2str(nanmean(val))])
        end
        plot([nanmean(val) nanmean(val)],[z(1) z(end)], '--','color',bluea,'linewidth',0.7)
        grid on
        ylim([0 42]); 
        if i>1
            xlim(xxplot(j,:));
        end
        if i==nm
            xlabel(uom{j})
        end
        if j==1
           ylabel({models{i};'Depth [m]'});
        else
           set(gca,'YTickLabel','')
        end
        set(gca,'YTick',[0:10:40])
    end
end
set(ha(10:27),'XTickLabel',''); %set(ha,'YTickLabel','')



tab2=readtable('tot_qual.csv');
c=0;
for i=1:nm
    for j=1:nmet
        c=c+1;
        idx=find(strcmp(table2array(tab2(:,1)),metrics{j}));
        val=tab2.(models{i})(idx);
        axes(ha(c));
        plot([nanmean(val) nanmean(val)],[z(1) z(end)], '--','color',reda,'linewidth',0.7)
    end
end
legend({'Depth specific','Depth averaged','Whole profile'})
print(gcf,'plot_metrics.png','-dpng','-r300')
print(gcf,'plot_metrics.pdf','-dpdf','-r300')


%%
figure('Units','centimeters','PaperSize',[15,12],'PaperPosition',[0 0 15 12])
ha = tight_subplot((nm-1),nmet,[.01 .02],[.12 .01],[.095 .01]);
c=0;
for i=2:nm
    for j=1:nmet
        c=c+1;
        idx=find(strcmp(tab.Model,models{i}) & strcmp(tab.Metric,metrics{j}));
        z=tab.Depth(idx);
        val=tab.Value(idx);
        axes(ha(c));
        plot(val,z,'-','color',blue,'linewidth',1.1); set(gca,'ydir','reverse'); hold on
        ii=find(val~=-Inf); val=val(ii);
%         plot([nanmean(val) nanmean(val)],[z(1) z(end)], '--','color',bluea,'linewidth',0.7)
        grid on
        ylim([0 42]); 
        if i>1
            xlim(xxplot(j,:));
        end
        if i==nm
            xlabel(uom{j})
        end
        if j==1
           ylabel({models{i};'Depth [m]'});
        else
           set(gca,'YTickLabel','')
        end
        set(gca,'YTick',[0:10:40])
    end
end
set(ha(1:18),'XTickLabel',''); %set(ha,'YTickLabel','')

c=0;
for i=2:nm
    for j=1:nmet
        c=c+1;
        idx=find(strcmp(table2array(tab2(:,1)),metrics{j}));
        val=tab2.(models{i})(idx);
        axes(ha(c));
        plot([nanmean(val) nanmean(val)],[z(1) z(end)], '-','color',red,'linewidth',1.1)
    end
end
legend({'Depth specific','Whole profile'})
print(gcf,'plot_metrics_def.png','-dpng','-r300')
print(gcf,'plot_metrics_def.pdf','-dpdf','-r300')

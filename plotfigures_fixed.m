function num_plot=plotfigures_fixed(conditionx,modex,n_segment,k,j0,indxx,plto3d,pltr3dr,pltr3d,pltperror2d,pltperror3dx,pltrecon_error3d,std1,dbfolders)
% plot the figures for fixed gantry therapy.
ii=0;
if conditionx
    ii=ii+1;
    if modex==1
        figure(n_segment), subplot(4,1,1)
        set(gcf,'Position',[5,30,1400,800])
        plot(1:length(indxx),plto3d(1:length(indxx),1),'r',... 
            1:length(indxx),pltr3d(1:length(indxx),1),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'), xlim([1 length(indxx)]),...
            legend('original SI','estimated SI'); ...
            title(strcat('original, estimated, and predicted positions of - (j0)',...
            dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)
        
        figure(n_segment), subplot(4,1,2)
        plot(1:length(indxx),plto3d(1:length(indxx),2),'r',... 
            1:length(indxx),pltr3d(1:length(indxx),2),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'), xlim([1 length(indxx)]),...
            legend('original LR','estimated LR'); 
        
        figure(n_segment), subplot(4,1,3)
        plot(1:length(indxx),plto3d(1:length(indxx),3),'r',... 
            1:length(indxx),pltr3d(1:length(indxx),3),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'),xlim([1 length(indxx)]),...
            legend('original AP','estimated AP'); 
        figure(n_segment), subplot(4,1,4)
        plot(1:length(indxx),pltrecon_error3d(1:length(indxx)),'m','LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Error(mm)'), xlabel('time in second'),...
            xlim([1 length(indxx)]),legend('3D estimation error');    
    elseif modex==2
        figure(n_segment), subplot(4,1,1)
        set(gcf,'Position',[5,30,1400,800])
        plot(1:length(indxx),plto3d(1:length(indxx),1),'r',... 
            1:length(indxx),pltr3dr(1:length(indxx),1),'b', 'LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'), xlim([1 length(indxx)]),...
            legend('original SI','predicted SI'); ...
            title(strcat('original, estimated, and predicted positions of - (j0)',...
            dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)        
        figure(n_segment), subplot(4,1,2)
        plot(1:length(indxx),plto3d(1:length(indxx),2),'r',... 
            1:length(indxx),pltr3dr(1:length(indxx),2),'b',...
            'LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'), xlim([1 length(indxx)]),...
            legend('original LR','predicted LR');         
        figure(n_segment), subplot(4,1,3)
        plot(1:length(indxx),plto3d(1:length(indxx),3),'r',... 
            1:length(indxx),pltr3dr(1:length(indxx),3),'b','LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Position(mm)'),xlim([1 length(indxx)]),...
            legend('original AP','predicted AP'); 
        figure(n_segment), subplot(4,1,4)            
        plot(1:length(indxx),pltperror2d(1:length(indxx)),'m',...
            1:length(indxx),pltperror3dx(1:length(indxx)),'k',... 
            'LineWidth',1.5)
        set(gca,'FontSize',16)
        set(gca,'XTick',[131 173 303 345 475 517 647 689 819],'XTickLabel',[20;70;90;140;160;210;230;280;300])
        grid()
        ylabel('Error(mm)'), xlabel('time in second'),...
            xlim([1 length(indxx)]),legend('2D combined error(2D beam¡¯s eye view error)','3D combined error'); 
    end
end
num_plot=ii;
return        

        
        
        
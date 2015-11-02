function plotfigures(i0,conditionx,modex,is_block,n_segment,points_period,xcoord,so3d,r3dr,r3d,perror2d,perror3dx,recon_error3d,runlength,ang_incrt,centerx0,dir_vector)
% to plot figure according to the condition and modex
% return the number of figures that satisfy the condition
ii=0;
j0=1;
if conditionx
    ii=ii+1;
    if modex==1
        figure(n_segment), subplot(4,1,1)
        set(gcf,'Position',[5,30,1400,800])
        plot(xcoord(i0:points_period),so3d(i0:points_period,1),'r',...  
            xcoord(i0:points_period),r3d(i0:points_period,1),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original SI', 'estimated SI');
%             title(strcat('original and estimated positions of - (j0)',...
%             dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)
        
        figure(n_segment), subplot(4,1,2)
        plot(xcoord(i0:points_period),so3d(i0:points_period,2),'r',...       
            xcoord(i0:points_period),r3d(i0:points_period,2),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original LR', 'estimated LR'); 
        
        figure(n_segment), subplot(4,1,3)
        plot(xcoord(i0:points_period),so3d(i0:points_period,3),'r',... 
            xcoord(i0:points_period),r3d(i0:points_period,3),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'),xlim([1 359]),...
            legend('original AP','estimated AP'); 
        if is_block==0
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),recon_error3d(i0:points_period),'m','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),...
                xlim([1 359]),legend('3D estimation error');   
        else
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),recon_error3d(i0:points_period),'m',...
                xcoord(i0:points_period),is_kv(i0:points_period,1),'c','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),... 
                xlim([1 359]),legend('3D estimation error','kV imaging on');   
        end
        figure(n_segment+190);
        for nn=1:floor(points_period/runlength)-1
%         subplot(2,2,rem(nn,4)+1);
%         plot3(...r3d(i0:points_period,1),r3d(i0:points_period,2),r3d(i0:points_period,3),'b',...
%             so3d(i0+(nn-1)*runlength:i0+nn*runlength,1),so3d(i0+(nn-1)*runlength:i0+nn*runlength,2),so3d(i0+(nn-1)*runlength:i0+nn*runlength,3),'LineWidth',1.5)
        plot3(r3d(i0+(nn-1)*runlength:i0+nn*runlength,1),r3d(i0+(nn-1)*runlength:i0+nn*runlength,2),r3d(i0+(nn-1)*runlength:i0+nn*runlength,3),'b',...
            so3d(i0+(nn-1)*runlength:i0+nn*runlength,1),so3d(i0+(nn-1)*runlength:i0+nn*runlength,2),so3d(i0+(nn-1)*runlength:i0+nn*runlength,3),'k','LineWidth',1.5)
        hold on;
%         t=-5:0.5:5;x=dir_vector(1)*t+centerx0(1);y=dir_vector(2)*t+centerx0(2);z=dir_vector(3)*t+centerx0(3);plot3(x,y,z);
         t=-5:0.5:5;x=dir_vector(i0+(nn-1)*runlength,1)*t+centerx0(i0+(nn-1)*runlength,1);y=dir_vector(i0+(nn-1)*runlength,2)*t+centerx0(i0+(nn-1)*runlength,2);z=dir_vector(i0+(nn-1)*runlength,3)*t+centerx0(i0+(nn-1)*runlength,3);plot3(x,y,z,'g','LineWidth',1.5);
         plot3(r3d(1:runlength,1),r3d(1:runlength,2),r3d(1:runlength,3),'r');
                set(gca,'FontSize',16)
        xlabel('SI position(mm)'),ylabel('LR position(mm)'),zlabel('AP position(mm)'),title((i0+(nn-1)*runlength)*ang_incrt);
%         legend('estimated','original');
        grid on
        pause
        cla(gca)
        end
%         title(strcat('original and estimated positions of - (j0)',...
%             dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)

    elseif modex==2
        figure(n_segment), subplot(4,1,1)
        set(gcf,'Position',[5,30,1400,800])
        plot(xcoord(i0:points_period),so3d(i0:points_period,1),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,1),'b','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original SI','predicted SI'); 
%             title(strcat('original and predicted positions of - (j0)',...
%             dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)
        
        figure(n_segment), subplot(4,1,2)
        plot(xcoord(i0:points_period),so3d(i0:points_period,2),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,2),'b','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original LR','predicted LR');         
        figure(n_segment), subplot(4,1,3)
        plot(xcoord(i0:points_period),so3d(i0:points_period,3),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,3),'b','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'),xlim([1 359]),...
            legend('original AP','predicted AP'); 
        if is_block==0
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),perror2d(i0:points_period),'m',...
                xcoord(i0:points_period),perror3dx(i0:points_period),'k','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),...
                xlim([1 359]),legend('2D combined error(2D beam¡¯s eye view error)','3D combined error');   
        else
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),perror2d(i0:points_period),'m',...
                xcoord(i0:points_period),perror3dx(i0:points_period),'k',... 
                xcoord(i0:points_period),is_kv(i0:points_period,1),'c','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),...
                xlim([1 359]),legend('2D combined error(2D beam¡¯s eye view error)','3D combined error','kV imaging on');   
        end
    elseif modex==3
        figure(n_segment), subplot(4,1,1)
        set(gcf,'Position',[5,30,1400,800])
        plot(xcoord(i0:points_period),so3d(i0:points_period,1),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,1),'b',...
            xcoord(i0:points_period),r3d(i0:points_period,1),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original SI','predicted SI', 'estimated SI');
%             title(strcat('original, estimated, and predicted positions of - (j0)',...
%             dbfolders(k).name,', j0=',int2str(j0), ', 2D rms error=', num2str(std1)),'FontSize',18)
        
        figure(n_segment), subplot(4,1,2)
        plot(xcoord(i0:points_period),so3d(i0:points_period,2),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,2),'b',...
            xcoord(i0:points_period),r3d(i0:points_period,2),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'), xlim([1 359]),...
            legend('original LR','predicted LR', 'estimated LR');
        
        figure(n_segment), subplot(4,1,3)
        plot(xcoord(i0:points_period),so3d(i0:points_period,3),'r',... 
            xcoord(i0:points_period),r3dr(i0:points_period,3),'b',...
            xcoord(i0:points_period),r3d(i0:points_period,3),'g','LineWidth',1.5)
        set(gca,'FontSize',16)
        grid()
        ylabel('Position(mm)'),xlim([1 359]),...
            legend('original AP','predicted AP', 'estimated AP');
        if is_block==0
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),perror2d(i0:points_period),'m',...
                xcoord(i0:points_period),perror3dx(i0:points_period),'k',...
                xcoord(i0:points_period),recon_error3d(i0:points_period),'c','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),...
                xlim([1 359]),legend('2D combined error(2D beam¡¯s eye view error)','3D combined error','3D estimation error');   
        else
            figure(n_segment), subplot(4,1,4)
            plot(xcoord(i0:points_period),perror2d(i0:points_period),'m',...
                xcoord(i0:points_period),perror3dx(i0:points_period),'k',... 
                xcoord(i0:points_period),recon_error3d(i0:points_period),'c',...
                xcoord(i0:points_period),is_kv(i0:points_period,1),'y','LineWidth',1.5)
            set(gca,'FontSize',16)
            grid()
            ylabel('Error(mm)'), xlabel('angles in degree'),...
                xlim([1 359]),legend('2D combined error(2D beam¡¯s eye view error)','3D combined error','3D estimation error','kV imaging on');   
        end
    end % of "if modex==..."
end % of 
num_plot=ii;
return
    
function [r3dx,centerx,kx,ko,ym,linearity,isplane,Bc,cost]=proj_backproj_new(seg_err,...
    marginx,so3d,r3d,jx,gAng,runlength,sid,stdx)
% this function projects the actual points to the MV imager and 
% reconstruct them

% r3dx returns the new reconstruction result.
% is_discon signals the discontinuity in the original data.
% is_jump signals the large velocity in the actual 3D motion.

% seg_err is the detection error(=0.5mm)
% k is the no. of fraction.
% r3d is the estimated 3D positions
% jx is the no. of point to be estimated.
% delta_angle is the angle of gantry rotation during the first 4 seconds.

%% Initialize the parameters. 
% centerx=[0 0 0];
kx=[0 0 0];
ko=[0 0 0];
ym=[0 0];
isplane=0;
linearity=1;
xfit=zeros(jx,3);
r3d0=zeros(jx,3);

% it has been compared to 1:jx-1. The following is better.
% for ix=1:runlength     
% %     r3d0(ix,1)=r3d(ix,2);
% %     r3d0(ix,2)=r3d(ix,1);
% %     r3d0(ix,3)=stdx-r3d(ix,3);
% end

p2d=[0 0 0];
wc3=[0 0 0];
prev = [0 0 0];
% if the point data in the database is continuous(based on the time 
% information), perform estimation, flag it otherwise.
o_3d=[0 0 0];
o_3d(:)=so3d(jx,:);

m=find(r3d(1:runlength,4)==1);
n=find(r3d(1:runlength,4)==0);
r3d11=r3d(r3d(1:runlength,4)==1,1:3);
r3d22=r3d(r3d(1:runlength,4)==0,1:3);
[coeff0,score0,~] = princomp(r3d(1:runlength,1:3));
% centerx0 = mean(r3d(1:runlength,1:3),1);
% xfit0=repmat(centerx0,runlength,1) + score0(:,1)*coeff0(:,1)';
% ymin0=min(xfit0(:,1));
% ymax0=max(xfit0(:,1));
% ymin0=min(score0(:,1));
% ymax0=max(score0(:,1));
for ix=1:runlength
   phase(ix)=(xfit(ix,1)-ymin)/(ymax-ymin);
%    phase(ix)=(score0(ix,1)-ymin0)/(ymax0-ymin0);
end
phase1=phase(m);
phase2=phase(n);   
% for ix=1:size(r3d11,1)
% curv1(ix)=curvature(r3d11,ix);
% end
% for ix=1:size(r3d22,1)
% curv2(ix)=curvature(r3d22,ix);
% end

%% project forward
r3d1=r3d(r3d(1:jx-1,4)==1,1:3);
r3d2=r3d(r3d(1:jx-1,4)==0,1:3);

p2d(1:3)=proj([o_3d(1) o_3d(2) o_3d(3)],gAng,sid,stdx);
[coeff,score,roots] = princomp(r3d(1:jx-1,1:3)) ;% only earlier points can be referred to in the calculation
r3d
[~,~,roots1] = princomp(r3d1(1:end,1:3)); % only earlier points can be referred to in the calculation
[~,~,roots2] = princomp(r3d2(1:end,1:3)); % only earlier points can be referred to in the calculation
centerx = [mean(r3d(1:jx-1,1:3),1);mean(r3d1(1:end,1:3),1);mean(r3d2(1:end,1:3),1)];
xfit=repmat(centerx(1,:),jx-1,1) + score(:,1)*coeff(:,1)';
% figure(100);
% plot3(xfit(1:jx-1,1),xfit(1:jx-1,2),xfit(1:jx-1,3))
% t=-5:0.5:5;x=coeff(1,1)*t+centerx(1,1);y=coeff(2,1)*t+centerx(1,2);z=coeff(3,1)*t+centerx(1,3);
% hold on
% plot3(x,y,z,'g','LineWidth',1.5);
% pause
linearity1 = roots1(1)/sum(roots1);
linearity2 = roots2(1)/sum(roots2);

% if linearity>0.001   % a parameter can be adjusted, it is no use though.
%     isplane=0;
% else
%     isplane=1;
% end
dir_vector = coeff(:,1)';
% diro_vector = coeffo(:,1)';
kx(:)=dir_vector(:);
% ko(:)=diro_vector(:);
%% add segmentation error to the data
p2d(1)=p2d(1) + seg_err * randn();
p2d(2)=p2d(2) + seg_err * randn();
ymin=min(xfit(:,1));
ymax=max(xfit(:,1));  
ym=[ymin ymax];
Bc=[0 0 0];
marginx = marginx*(ymax-ymin);

prev = r3d(jx-1,1:3) - r3d(jx-2,1:3);
prep = r3d(jx-1,1:3);
%% backproject according to whether it's a plane or not.
% if isplane==0
    [wc3,B,cost]=back_proj_line(centerx(1,:),kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
    r3dx(1)= wc3(2)
    r3dx(2)= wc3(1)
    r3dx(3)= -wc3(3)+stdx%将坐标原点移至旋转框架的中心，即世界坐标系的原点
    current_position=bsxfun(@minus,r3dx(1:3),mean([r3d(1:jx-1,1:3);r3dx(1:3)],1))*coeff(:,1)/norm(bsxfun(@minus,r3dx(1:3),mean([r3d(1:jx-1,1:3);r3dx(1:3)],1)))
    pre_position=(xfit(jx-1,:)-centerx(1,:))*coeff(:,1)
    if (current_position-pre_position)>=0.001   
%         [wc3 B cost]=back_proj_line(centerx(2,:),kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
%         r3d(jx,1)= wc3(2);
%         r3d(jx,2)= wc3(1);
%         r3d(jx,3)= -wc3(3)+stdx;%将坐标原点移至旋转框架的中心，即世界坐标系的原点 
%         current_score=bsxfun(@minus,r3dx,mean(r3d(1:jx-1,1:3),1))*coeff(:,1)
        if current_position>=ymin0 && current_position<=ymax0
           dif=abs(r3d1(:,1)-r3dx(1));
           dif_sorted=sort(dif)
           kxx=r3d1(dif==dif_sorted(1),1:3)-r3d1(dif==dif_sorted(2),1:3);
           [wc3 B cost]=back_proj_line(r3d11(dif==dif_sorted(1),1:3),kxx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
        end    
        r3dx(4)= 1;
        
%         closest_curvature=curv1(dif==min(dif))
%         [wc3,lambda]=brent(r3d,jx,B,closest_curvature,stdx)
%          centerx=centerx(2,:);
    elseif (current_position-pre_position)<=0.001
%         [wc3 B cost]=back_proj_line(centerx(3,:),kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
%         r3d(jx,1)= wc3(2);
%         r3d(jx,2)= wc3(1);
%         r3d(jx,3)= -wc3(3)+stdx;%将坐标原点移至旋转框架的中心，即世界坐标系的原点
%         current_score=bsxfun(@minus,r3dx,mean(r3d(1:jx-1,1:3),1))*coeff(:,1)
        if current_position>=ymin0 && current_position<=ymax0
           dif=abs(phase2-(current_position-ymin0)/(ymax0-ymin0));
           dif_sorted=sort(dif)
           kxx=r3d22(dif==dif_sorted(1),1:3)-r3d22(dif==dif_sorted(2),1:3)
           [wc3 B cost]=back_proj_line(r3d22(dif==dif_sorted(1),1:3),kxx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
        end
        r3dx(4)=0;
%         closest_curvature=curv2(dif==min(dif))
%         [wc3,lambda]=brent(r3d,jx,B,closest_curvature,stdx)
%          centerx=centerx(3,:);
    else
        if r3d(jx-1,4)==1
            if current_position>=ymin0 && current_position<=ymax0
            dif=abs(phase1-(current_position-ymin0)/(ymax0-ymin0))
            dif_sorted=sort(dif)
            kx=r3d11(dif==dif_sorted(1),1:3)-r3d11(dif==dif_sorted(2),1:3);
            [wc3 B cost]=back_proj_line(r3d11(dif==dif_sorted(1),1:3),kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
            end
            r3dx(4)= 1;
        else
            if current_position>=ymin0 && current_position<=ymax0
            dif=abs(phase2-(current_position-ymin0)/(ymax0-ymin0));
            dif_sorted=sort(dif)
            kx=r3d22(dif==dif_sorted(1),1:3)-r3d22(dif==dif_sorted(2),1:3)
            [wc3 B cost]=back_proj_line(r3d22(dif==dif_sorted(1),1:3),kx,prev,prep,marginx,ymin,ymax,stdx,sid,gAng,p2d(1),p2d(2));
            end
            r3dx(4)= 0;
        end
    end
    r3dx(1)= wc3(2);
    r3dx(2)= wc3(1);
    r3dx(3)= -wc3(3)+stdx;%将坐标原点移至旋转框架的中心，即世界坐标系的原点
    Bc(1) = B(2);
    Bc(2) = B(1);
    Bc(3) = -B(3);
    Bc = Bc/norm(Bc);
    centerx=centerx(1,:);
% %     pause
% else
%     wc3=back_proj_plane(centerx,kx,r3d0,stdx,sid,gAng,p2d(1),p2d(2));
%     r3dx(1)= wc3(2);
%     r3dx(2)= wc3(1);
%     r3dx(3)= -wc3(3)+stdx;
% end
% xx=range(so3d(1:jx,1:3),1);
% range_o3d=xx(1,:); % ??
return
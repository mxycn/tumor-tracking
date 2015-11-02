function [s cost]=find_s(O3,B3,c,kp,prev,prep,marginx,y_min,y_max,r,factorv)
% Use cost function and iterative method to find the most suitable s, hence the best point.
% Express O3 and B3 in the common coordinate system
O(1) = O3(2);
O(2) = O3(1);
O(3) =r- O3(3);
B(1) = B3(2);
B(2) = B3(1);
B(3) = - B3(3);
cos2 = B*kp'/norm(B);
cos2 = cos2*cos2;
    % Get close to the line
    % cost function  = 
    %            (1-cos2)*distance to line^2 +
    %            factorv*cos2*(currentpoint-prep-prev)^2
    % distance between a point and a line, refer to 
    % Handbook of mathematics. P 218.
%     distance2 = ...
%     1/norm(kp)^2*(((O(1)+s*B(1)-c(1))*kp(2)-(O(2)+s*B(2)-c(2))*kp(1))^2 + ...
%     ((O(2)+s*B(2)-c(2))*kp(3)-(O(3)+s*B(3)-c(3))*kp(2))^2 + ...
%     ((O(3)+s*B(3)-c(3))*kp(1)-(O(1)+s*B(1)-c(1))*kp(3))^2);
    
    % (currentpoint-prep-prev)^2 = 
    % norm(O+s*B-prep-prev)^2
%     term2 = factorv*norm(O+s*B-prep-prev)^2;
    
    % cost = s^2*((B(1)*kp(2)-B(2)*kp(1))^2+(B(2)*kp(3)-B(3)*kp(2))^2+...
    %                    (B(3)*kp(1)-B(1)*kp(3))^2)/(kp(1)^2+kp(2)^2+kp(3)^2)  + ...
    %         s*((B(1)*kp(2)-B(2)*kp(1))*((O(1)-c(1))*kp(2)-(O(2)-c(2))*kp(1)) + ...
    %            (B(2)*kp(3)-B(3)*kp(2))*((O(2)-c(2))*kp(3)-(O(3)-c(3))*kp(2))  + ...
    %            (B(3)*kp(1)-B(1)*kp(3))*((O(3)-c(3))*kp(1)-(O(1)-c(1))*kp(3)))/(kp(1)^2+kp(2)^2+kp(3)^2) 
    %            + ... + ...
    %         s^2*(B(1)^2+B(2)^2+B(3)^2)*cos2*factorv + ...
    %         s*B*(O-prep-prev)' + ...

%     aa = ((B(1)*kp(2)-B(2)*kp(1))^2+(B(2)*kp(3)-B(3)*kp(2))^2+...
%         (B(3)*kp(1)-B(1)*kp(3))^2) + ...
%         (B(1)^2+B(2)^2+B(3)^2)*factorv*cos2;
    aa1 = (B(1)*kp(2)-B(2)*kp(1))^2;
    aa2 = (B(2)*kp(3)-B(3)*kp(2))^2;
    aa3 = (B(3)*kp(1)-B(1)*kp(3))^2;
    aa4 = (B(1)^2+B(2)^2+B(3)^2)*factorv*cos2;
    aa = aa1 + aa2 + aa3+ aa4;    
%     bb = ((B(1)*kp(2)-B(2)*kp(1))*((O(1)-c(1))*kp(2)-(O(2)-c(2))*kp(1)) + ...
%         (B(2)*kp(3)-B(3)*kp(2))*((O(2)-c(2))*kp(3)-(O(3)-c(3))*kp(2))  + ...
%         (B(3)*kp(1)-B(1)*kp(3))*((O(3)-c(3))*kp(1)-(O(1)-c(1))*kp(3))) + ...
%         B*(O-prep-prev)'*factorv*cos2;
    bb1 = (B(1)*kp(2)-B(2)*kp(1))*((O(1)-c(1))*kp(2)-(O(2)-c(2))*kp(1));
    bb2 = (B(2)*kp(3)-B(3)*kp(2))*((O(2)-c(2))*kp(3)-(O(3)-c(3))*kp(2));
    bb3 = (B(3)*kp(1)-B(1)*kp(3))*((O(3)-c(3))*kp(1)-(O(1)-c(1))*kp(3));
    bb4 = B*(O-prep-prev)'*factorv*cos2;
    bb = bb1 + bb2 + bb3 + bb4;
     s = -bb/aa;
%     distance2 = ...
%     1/norm(kp)^2*(((O(1)+s*B(1)-c(1))*kp(2)-(O(2)+s*B(2)-c(2))*kp(1))^2 + ...
%     ((O(2)+s*B(2)-c(2))*kp(3)-(O(3)+s*B(3)-c(3))*kp(2))^2 + ...
%     ((O(3)+s*B(3)-c(3))*kp(1)-(O(1)+s*B(1)-c(1))*kp(3))^2);
    distance21 = ((O(1)+s*B(1)-c(1))*kp(2)-(O(2)+s*B(2)-c(2))*kp(1))^2;
    distance22 = ((O(2)+s*B(2)-c(2))*kp(3)-(O(3)+s*B(3)-c(3))*kp(2))^2;
    distance23 = ((O(3)+s*B(3)-c(3))*kp(1)-(O(1)+s*B(1)-c(1))*kp(3))^2;
    distance2 = (distance21+ distance22+ distance23)/norm(kp)^2;
    term2 = factorv*norm(O+s*B-prep-prev)^2;
%     Get close to the nearest point
    cost = (1-cos2)*distance2 + term2*cos2;
    P = O3+s*B3;
        % if the point is close to the check line extension, instead of the
    % checkline segment, then use another method.
    
    k1(1)=kp(2)/kp(1);
    k1(2)=-kp(2)/kp(1)*c(1)+c(2);
    k2(1)=-kp(3)/kp(1);
    k2(2)=kp(3)/kp(1)*c(1)+r-c(3);
    bbx=P(2)-k1(1)*(k1(2)-P(1))-k2(1)*(k2(2)-P(3));
    aax=k1(1)*k1(1)+k2(1)*k2(1)+1;    
    y=bbx/aax; % the y-coordinate of the point on the trajectory line that is closest to the target 
 % negative means conservative, positive means following trend     
     if y_min-y>marginx*abs(kp(1))             
        xmin=k1(1)*y_min+k1(2);
        ymin=y_min;
        zmin=k2(1)*y_min+k2(2);
        a1=O3(1)-xmin;
        a2=O3(2)-ymin;
        a3=O3(3)-zmin;
        b1=B3(1);
        b2=B3(2);
        b3=B3(3);
        s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
    elseif y-y_max>marginx*abs(kp(1))
%     elseif y-y_max>margin 
        xmax=k1(1)*y_max+k1(2);
        ymax=y_max;
        zmax=k2(1)*y_max+k2(2);
        a1=O3(1)-xmax;
        a2=O3(2)-ymax;
        a3=O3(3)-zmax;
        b1=B3(1);
        b2=B3(2);
        b3=B3(3);
        s=-(a1*b1+a2*b2+a3*b3)/(b1*b1+b2*b2+b3*b3);
     end   
    
    
    % distance2 =
    % 1/norm(kp)^2*(((O(1)+s*B(1)-c(1))*kp(2)-(O(2)+s*B(2)-c(2))*kp(1))^2 +
    % ((O(2)+s*B(2)-c(2))*kp(3)-(O(3)+s*B(3)-c(3))*kp(2))^2 +
    % ((O(3)+s*B(3)-c(3))*kp(1)-(O(1)+s*B(1)-c(1))*kp(3))^2)
    
    % (currentpoint-prep-prev)^2 = 
    % norm(O+s*B-prep-prev)^2
    
    % cost = s^2*((B(1)*kp(2)-B(2)*kp(1))^2+(B(2)*kp(3)-B3)*kp(2))^2+...
    %                    (B(3)*kp(1)-B(1)*kp(3))^2)/(kp(1)^2+kp(2)^2+kp(3)^2)  + ...
    %         s*((B(1)*kp(2)-B(2)*kp(1))*((O(1)-c(1))*kp(2)-(O(2)-c(2))*kp(1)) + ...
    %            (B(2)*kp(3)-B(3)*kp(2))*((O(2)-c(2))*kp(3)-(O(3)-c(3))*kp(2))  + ...
    %            (B(3)*kp(1)-B(1)*kp(3))*((O(3)-c(3))*kp(1)-(O(1)-c(1))*kp(3)))/(kp(1)^2+kp(2)^2+kp(3)^2) 
    %            + ... + ...
    %         s^2*(B(1)^2+B(2)^2+B(3)^2) + ...
    %         s*B*(O-prep-prev)' + ...
    
    
    
    

    
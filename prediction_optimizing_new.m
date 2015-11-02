function r3dr=prediction_optimizing_new(jx,r3dp,stdx,centerx,kx,linearity);
% use the linear correlation between the motions of different direction to
% optimize the prediction result
% input: jx-the no. of the current point; r3dp-the predicted position
% before optimization; stdx

kk=linearity*linearity*linearity;
k1=[0 0];
k2=[0 0];
% if abs(kx(2))<0.001
%     pause(0.1);
% end
k1(1)=kx(2)/kx(1);
k1(2)=-kx(2)/kx(1)*centerx(1)+centerx(2);%轨迹在XY平面投影与Y轴的交点的Y坐标
k2(1)=-kx(3)/kx(1);
k2(2)=kx(3)/kx(1)*centerx(1)+stdx-centerx(3);%轨迹在XZ平面投影与Z轴的交点的Z坐标，此时坐标原点为MV_sourece,Z轴正方向为射线方向

% kk=0;
r3dr=[r3dp(jx,1) r3dp(jx,2) r3dp(jx,3)];
aa=r3dp(jx,1)-k1(1)*(k1(2)-r3dp(jx,2))-k2(1)*(k2(2)-stdx+r3dp(jx,3));
bb=1+k1(1)*k1(1)+k2(1)*k2(1);
xx=aa/bb;
r3dr(1)=kk*xx+(1-kk)*r3dp(jx,1);
r3dr(2)=kk*(k1(1)*xx+k1(2))+(1-kk)*r3dp(jx,2);
r3dr(3)=kk*(stdx-(k2(1)*xx+k2(2)))+(1-kk)*r3dp(jx,3);
r3dr(4)=r3dp(jx,4);
return
% proj.m
% project the world coord in mm of a target to MV imager at any gantry angle
% world coord definition: X+ is inferior, Y+ is left, and Z+ is anterior
% 07/03/10 suppose fixed gantry angle

% revised by Huagang Yan, Aug 12,2010

function pc = proj(wc, gAng, F_mv, R_mv)

%%% INPUT: refer to "GetTRdoc.pdf", page 42.
% wc(3) - an array. 3D world coord of a target: X+ is inferior, Y+ is left, and Z+ is anterior
% gAng - gantry angle in units of deg the starting ray is in z+ direction.
%     F_mv = MV source-to-MV imager distance
%     R_mv = the projection of MV source-to-the origin distance on the ratational plane(YZ plane), 
%                 where the origin is in coincidence with the rotation center on the YZ plane.
%%% OUTPUT
% pc - pixel coord of the projection of target on the MV imager 2xN array, the center is the origin
wcc=wc+[0 0 -R_mv];

wcp=[0 0 0];
wcp(1)=wcc(1);
wcp(2)=wcc(2)*cos(gAng/180*pi)+(wcc(3)+R_mv)*sin(gAng/180*pi);
wcp(3)=-wcc(2)*sin(gAng/180*pi)+(wcc(3)+R_mv)*cos(gAng/180*pi)-R_mv;

ratiop=-F_mv/wcp(3);

pc = [0 0 0];
pc(1)= wcp(2)*ratiop;
pc(2)= wcp(1)*ratiop;
pc(3) = ratiop;
return



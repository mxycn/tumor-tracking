function [gAng,is_do]=assign_angle_fixed(jj,j0,ko3d,angleset,onofftime)
% assign gantry angle according to treatment protocol.
if ko3d(jj,4)-ko3d(j0,4)>=onofftime(1) & ko3d(jj,4)-ko3d(j0,4)<onofftime(2)
    gAng=angleset(1);
    is_do=1;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(2) & ko3d(jj,4)-ko3d(j0,4)<onofftime(3)
    gAng=NaN;
    is_do=0;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(3) & ko3d(jj,4)-ko3d(j0,4)<onofftime(4)
    gAng=angleset(2);
    is_do=1;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(4) & ko3d(jj,4)-ko3d(j0,4)<onofftime(5)
    gAng=NaN;
    is_do=0;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(5) & ko3d(jj,4)-ko3d(j0,4)<onofftime(6)
    gAng=angleset(3);
    is_do=1;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(6) & ko3d(jj,4)-ko3d(j0,4)<onofftime(7)
    gAng=NaN;
    is_do=0;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(7) & ko3d(jj,4)-ko3d(j0,4)<onofftime(8)
    gAng=angleset(4);
    is_do=1;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(8) & ko3d(jj,4)-ko3d(j0,4)<onofftime(9)
    gAng=NaN;
    is_do=0;
elseif ko3d(jj,4)-ko3d(j0,4)>=onofftime(9) & ko3d(jj,4)-ko3d(j0,4)<onofftime(10)
    gAng=angleset(5);
    is_do=1;
else
    gAng=NaN;
    is_do=0;
end
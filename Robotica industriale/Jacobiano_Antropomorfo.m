function [J]=Jacobiano_Antropomorfo(L,Q)
L1=L(1);
L2=L(2);
L3=L(3);
c1=cos(Q(1));c2=cos(Q(2));s1=sin(Q(1));s2=sin(Q(2));
c23=cos(Q(2)+Q(3));s23=sin(Q(2)+Q(3));

J=[-s1*(L2*c2+L3*c23),-c1*(L2*s2+L3*s23),-L3*s23*c1;
    c1*(L2*c2+L3*c23),-s1*(L2*s2+L3*s23),-L3*s23*s1;
    0                , L2*c2+L3*c23     , L3*c23;];
end
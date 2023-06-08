function [Q]=Antropomorfo_Cin_Inv(L,P)
%antropomorfo libro Siciliano
Px=P(1);
Py=P(2);
Pz=P(3);
%
%test spazio di lavoro
%
d=norm(P);
if d<=(L(2)+L(3)) 
q1=atan2(Py,Px);
c1=cos(q1);
s1=sin(q1);
%
c3=((Px)^2+Py^2+(Pz)^2-L(2)^2-L(3)^2)/(2*L(2)*L(3));
s3=sqrt(1-c3^2);
q3=atan2(s3,c3);
%compute Q(2)
A=[L(2)+L(3)*c3 L(3)*s3; -L(3)*s3 L(2)+L(3)*c3];
B=inv(A)*[Pz; sqrt(Px^2+Py^2)];
q2=atan2(B(1),B(2));

Q=[q1;q2;q3];
else
    Q=[0;0;0]
end
    

end
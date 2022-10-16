clc
clear all
close all
format short


syms q1 q2 q3 l0 l1 l2 l3 s1 s2 s3 s4 c1 c2 c3 c4 Px Py Pz
%%
rotz = [s1 c1 0 0; -c1 s1 0 0; 0 0 1 0; 0 0 0 1];
traz = [1 0 0 0; 0 1 0 0; 0 0 1 l0; 0 0 0 1];
trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
a01 = rotz*traz*trax*rotx
%%
rotz = [-s2 -c2 0 0; c2 -s2 0 0; 0 0 1 0; 0 0 0 1];
traz = [1 0 0 0; 0 1 0 0; 0 0 1 -l1; 0 0 0 1];
trax = [1 0 0 l2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
rotx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
a12 = rotz*traz*trax*rotx
a02 = a01*a12
%%
rotz = [c3  -s3 0 0; s3 c3 0 0; 0 0 1 0; 0 0 0 1];
traz = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
a23 = rotz*traz*trax*rotx
%%
rotz = [0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1];
traz = [1 0 0 0; 0 1 0 0; 0 0 1 l3; 0 0 0 1];
trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
rotx = [1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1];
a3e = rotz*traz*trax*rotx
%%
a = a01*a12*a23*a3e
%%
b0=[0;0;1;];
b1=a01(1:3,3);
b2=a02(1:3,3);
r0E=a(1:3,4);
r1E=a(1:3,4)-a01(1:3,4);
r2E=a(1:3,4)-a02(1:3,4);

JL1=cross(b0,r0E);
JL2=cross(b1,r1E);
JL3=cross(b2,r2E);

J=simplify([JL1,JL2,JL3; b0,b1,b2;])
JL=simplify(J(1:3,:))
det(JL)
JA=simplify(J(4:6,:))
det(JA)

inverted_J = simplify(inv(J(1:3,:)))

%%
eqns = [Px == cos(q1)*l1 + l3*(cos(q2)*cos(q3)*sin(q1) - sin(q1)*sin(q2)*sin(q3)) - l2*sin(q1)*sin(q2) Py == l1*sin(q1) - l3*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) + cos(q1)*l2*sin(q2) Pz == l0 + cos(q2)*l2 + l3*(cos(q2)*sin(q3) + cos(q3)*sin(q2))];

[s1,parameters,conditions] = solve(eqns,q1,"ReturnConditions",true)
[s2,parameters,conditions] = solve(eqns,q2,"ReturnConditions",true)
[s3,parameters,conditions] = solve(eqns,q3,"ReturnConditions",true)
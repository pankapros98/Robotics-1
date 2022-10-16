clc
clear all
close all

% *** Robot (kinematic) model parameters *** 
l0 = 10.0;  %% in cm 
l1 = 15.0;
l2 = 30.0;
l3 = l2;

%% *** sampling period *** 
% *** for the robot motion, kinematic simulation: 
dt = 0.001; %Tf = 0.001; i.e. 1 msec)   
%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
Tf=2.0; 	% 2sec Tf of motion 
t=dt:dt:2*Tf; 
%% initial & final end-point position & acceleration (cm,sec)
xA = 15;	
yA = -30;
zA = 40;

xB = 15;     
yB = -30;
zB = 50;

gx = 0.05;   
gy = 0.05;
gz = 0.05;

%% Desired Task-Space Trajectory %%%

xd  = double(trajectory(0,Tf,xA,xB,gx));
yd  = double(trajectory(0,Tf,yA,yB,gy));
zd  = double(trajectory(0,Tf,zA,zB,gz));
rxd = double(trajectory(Tf,2*Tf,xB,xA,-gx));
ryd = double(trajectory(Tf,2*Tf,yB,yA,-gy));
rzd = double(trajectory(Tf,2*Tf,zB,zA,-gz));

kmax = Tf/dt + 1;
x = zeros(length(t),1);
y= zeros(length(t),1);
z= zeros(length(t),1);
ux= zeros(length(t),1);
uy= zeros(length(t),1);
uz= zeros(length(t),1);

for k = 1:2*kmax-1
   time  = (k-1)*dt;
   if (time<Tf*0.1)
        x(k)  = polyval(xd(6:-1:1),time);
        ux(k) = polyval(polyder(xd(6:-1:1)),time);
        y(k)  = polyval(yd(6:-1:1),time);
        uy(k) = polyval(polyder(yd(6:-1:1)),time);
        z(k)  = polyval(zd(6:-1:1),time);
        uz(k) = polyval(polyder(zd(6:-1:1)),time);
   elseif (time<Tf*0.9)
        x(k)  = polyval(xd(8:-1:7),time);
        ux(k) = polyval(polyder(xd(8:-1:7)),time);
        y(k)  = polyval(yd(8:-1:7),time);
        uy(k) = polyval(polyder(yd(8:-1:7)),time);
        z(k)  = polyval(zd(8:-1:7),time);
        uz(k) = polyval(polyder(zd(8:-1:7)),time);
   elseif (time<Tf)
        x(k)  = polyval(xd(14:-1:9),time);
        ux(k) = polyval(polyder(xd(14:-1:9)),time);
        y(k)  = polyval(yd(14:-1:9),time);
        uy(k) = polyval(polyder(yd(14:-1:9)),time);
        z(k)  = polyval(zd(14:-1:9),time);
        uz(k) = polyval(polyder(zd(14:-1:9)),time);
   elseif (time<Tf*1.1)
        x(k)  = polyval(rxd(6:-1:1),time);
        ux(k) = polyval(polyder(rxd(6:-1:1)),time);
        y(k)  = polyval(ryd(6:-1:1),time);
        uy(k) = polyval(polyder(ryd(6:-1:1)),time);
        z(k)  = polyval(rzd(6:-1:1),time);
        uz(k) = polyval(polyder(rzd(6:-1:1)),time);
   elseif (time<Tf*1.9)
        x(k)  = polyval(rxd(8:-1:7),time);
        ux(k) = polyval(polyder(rxd(8:-1:7)),time);
        y(k)  = polyval(ryd(8:-1:7),time);
        uy(k) = polyval(polyder(ryd(8:-1:7)),time);
        z(k)  = polyval(rzd(8:-1:7),time);
        uz(k) = polyval(polyder(rzd(8:-1:7)),time);
   else
        x(k)  = polyval(rxd(14:-1:9),time);
        ux(k) = polyval(polyder(rxd(14:-1:9)),time);
        y(k)  = polyval(ryd(14:-1:9),time);
        uy(k) = polyval(polyder(ryd(14:-1:9)),time);
        z(k)  = polyval(rzd(14:-1:9),time);
        uz(k) = polyval(polyder(rzd(14:-1:9)),time);
   end
end
x  = x(1:length(t));   
y = y(1:length(t));  
z = z(1:length(t));
ux = ux(1:length(t)); 
uy = uy(1:length(t));
uz = uz(1:length(t));

% fig = figure;

subplot(2,3,1); 
plot(t,x); 
ylabel('x (cm)'); 
xlabel('time (sec)');
grid on

subplot(2,3,2); 
plot(t,y); 
ylabel('y (cm)'); 
xlabel('time (sec)');
grid on

subplot(2,3,3); 
plot(t,z); 
ylabel('z (cm)'); 
xlabel('time (sec)');
grid on

subplot(2,3,4); 
plot(t,ux); 
ylabel('ux (cm/sec)'); 
xlabel('time (sec)');
grid on

subplot(2,3,5); 
plot(t,uy); 
ylabel('uy (cm/sec)'); 
xlabel('time (sec)');
grid on

subplot(2,3,6); 
plot(t,uz); 
ylabel('uz (cm/sec)'); 
xlabel('time (sec)');
grid on

q3 = asin((x.^2 - y.^2 - l1^2 + (z - l0).^2 - l2^2 -l3^2)./(2* l2 * l3));
q2 = atan2(sqrt(x.^2 + y.^2 - l1^2), z-l0) + atan2(l3 .* cos(q3), l2 + l3.*cos(q3));
q1 = atan2(y, x) - atan2(l2.*sin(q2) - l3.*cos(q2).*cos(q3), l1);
qd1 = zeros(length(t),1); 
qd2= zeros(length(t),1); 
qd3= zeros(length(t),1);

for k = 1:length(t)
    j = inverse(l0,l1,l2,l3,q1(k),q2(k),q3(k));
    qd1(k) = j(1,1)*ux(k) + j(1,2)*uy(k);
    qd2(k) = j(2,1)*ux(k) + j(2,2)*uy(k);
    qd3(k) = j(3,1)*ux(k) + j(3,2)*uy(k);
end

fig2 = figure;

subplot(2,3,1); 
plot(t,q1);
ylabel('q1 (rad)'); 
xlabel('time (sec)');
grid on

subplot(2,3,2); 
plot(t,q2); 
ylabel('q2 (rad)'); 
xlabel('time (sec)');
grid on

subplot(2,3,3); 
plot(t,q3);
ylabel('q3 (rad)'); 
xlabel('time (sec)');
grid on

subplot(2,3,4); 
plot(t,qd1);
ylabel('qd1 (rad)'); 
xlabel('time (sec)');
grid on

subplot(2,3,5); 
plot(t,qd2); 
ylabel('qd2 (rad)'); 
xlabel('time (sec)');
grid on

subplot(2,3,6); 
plot(t,qd3);
ylabel('qd3 (rad)'); 
xlabel('time (sec)');
grid on

function X = trajectory(t0,tf,xA,xB,g)
    
    syms a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5
    t1 = t0+(tf-t0)*0.1; t2 = t0+(tf-t0)*0.9;
    eqn1 = a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5 - xA == 0;

    eqn2 = a1 + 2*a2*t0 + 3*a3*t0^2 + 4*a4*t0^3 + 5*a5*t0^4 == 0;

    eqn3 = 2*a2 + 6*a3*t0 + 12*a4*t0^2 + 20*a5*t0^3 - g == 0;

    eqn4 = a0 + a1*t1 + a2*t1^2 + a3*t1^3 + a4*t1^4 + a5*t1^5 - (b0 + b1*t1) == 0;

    eqn5 = a1 + 2*a2*t1 + 3*a3*t1^2 + 4*a4*t1^3 + 5*a5*t1^4 -b1 == 0;

    eqn6 = 2*a2 + 6*a3*t1 + 12*a4*t1^2 + 20*a5*t1^3 == 0;

    eqn7 = c0 + c1*t2 + c2*t2^2 + c3*t2^3 + c4*t2^4 + c5*t2^5 - (b0 + b1*t2) == 0;

    eqn8 = c1 + 2*c2*t2 + 3*c3*t2^2 + 4*c4*t2^3 + 5*c5*t2^4 - b1 == 0;
    
    eqn9 = 2*c2 + 6*c3*t2 + 12*c4*t2^2 + 20*c5*t2^3 == 0;

    eqn10 = c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5 - xB == 0;

    eqn11 = c1 + 2*c2*tf + 3*c3*tf^2 + 4*c4*tf^3 + 5*c5*tf^4 == 0;

    eqn12 = 2*c2 + 6*c3*tf + 12*c4*tf^2 + 20*c5*tf^3 +g == 0;

    eqn13 = b0+b1*t1 - (a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4 + a5*t0^5) - (0.5*(t1-t0)*b1) == 0;

    eqn14 = (c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4 + c5*tf^5) - (b0+b1*t2) - (0.5*(t1-t0)*b1) == 0;
    
    [A,B] = equationsToMatrix([ eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10 eqn11 eqn12 eqn13 eqn14 ] , [ a0 a1 a2 a3 a4 a5 b0 b1 c0 c1 c2 c3 c4 c5 ]);
    X = vpa(linsolve(A,B));
    
end

function Jinv = inverse(l0,l1,l2,l3,q1,q2,q3)

    rotz = [sin(q1) cos(q1) 0 0; -cos(q1) sin(q1) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 l0; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
    a01 = rotz*traz*trax*rotx;
%%
    rotz = [-sin(q2) -cos(q2) 0 0; cos(q2) -sin(q2) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 -l1; 0 0 0 1];
    trax = [1 0 0 l2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    a12 = rotz*traz*trax*rotx;
    a02 = a01*a12;
%%
    rotz = [cos(q3)  -sin(q3) 0 0; sin(q3) cos(q3) 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1];
    a23 = rotz*traz*trax*rotx;
%%
    rotz = [0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1];
    traz = [1 0 0 0; 0 1 0 0; 0 0 1 l3; 0 0 0 1];
    trax = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    rotx = [1 0 0 0; 0 0 1 0; 0 -1 0 0; 0 0 0 1];
    a3e = rotz*traz*trax*rotx;
%%
    a = a01*a12*a23*a3e;
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

    J= [JL1,JL2,JL3; b0,b1,b2];
    Jinv = inv(J(1:3,:));
end

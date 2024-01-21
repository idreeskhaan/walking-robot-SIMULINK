clc;
close all;
format short;

% 
m0= 0.036;
m1= 0.544;
m2= 0.302;
m3= 3.0;   %Torso
m4= 0.302; %upper leg 0.302
m5= 0.544; %lower leg 0.544
m6= 0.036; %Foot      0.036
SW_z = 0.05; %5 cm 
SH_y= 0.02;  %2cm
yH= 0.40;  %44cm max
T1 = 1;    %half step time
l0=0.045;
l1=0.180;
l2=0.220;
l4= 0.22; %upper leg 0.22
l5= 0.18; %lower leg 0.18
l6=0.045;  %Foot      0.045
g= 9.81;
lamda= sqrt(g/yH);
res=200;  %Resolution
y3=yH;

%define End effector orientation
t= linspace(-T1, T1, 5);
qE= [degtorad(-90) degtorad(-95) degtorad(-90) degtorad(-80) degtorad(-90)];
qE_c= polyfit(t, qE, 4);   %4th order fit

%Define Foot Trajectory
t=linspace(-T1, T1, res);
y6= SH_y*(1 - (t.^2)./(T1.^2));
z6= SW_z/2*((3.*t./T1) - (t.^3)./(T1.^3));


y6_c= polyfit(t, y6, 4);
z6_c= polyfit(t, z6, 4);

y6_dot_c = polyder(y6_c); 
z6_dot_c = polyder(z6_c); 

y6_2dot_c = polyder(y6_dot_c);
z6_2dot_c = polyder(z6_dot_c);

m3_t= m3+m2+m1;  %Rest of the mass of the robot
%Excitation function f(t) for TMIPM derived on the basis of ZMP Equation
m6_t= m4+m5+m6;
f_c= m6_t./(m3_t.*yH).*(g.*[0 0 0 0 z6_c] + conv(z6_c, [0 0 y6_2dot_c]) - conv(y6_c, [0 0 z6_2dot_c]));

%solve the ODE using bvp4c to get torso trajectory z3

t=linspace(-T1, T1, res);
solinit = bvpinit(t, @guess);
sol = bvp4c(@(t,y) [y(2), polyval(f_c,t) + lamda.^2.*y(1)], @(ya, yb) [ya(1) + SW_z/2, yb(1) - SW_z/2], solinit);
z3= sol.y(1,:);

% Convert z3 into a polynomial
t=linspace(-T1, T1, res);
z3_c= polyfit(t, z3, 5);
figure;
plot(t, polyval(z3_c,t), 'color', 'b', 'LineWidth', 2);
xlabel('t');
ylabel('z3 in m (TMIPM)');


% Inverse Kinematics
k=1;
while k>=1

j=1;    
for t=linspace(-T1, T1, res) 
    
%End effector location relative to Torso in World Frame
%px(j)= SW_x/2*((3.*t./T1) - (t.^3)./(T1.^3)); 
px(j)=0;
py(j)=polyval(y6_c,t)-yH;
pz(j)=polyval(z6_c,t)-polyval(z3_c,t);
%pz(j)=0;
phi(j) = polyval(qE_c,t); 

R= RY(-pi/2)*RZ(phi(j));
r11(j)= R(1,1); r12(j)= R(1,2); r13(j)= R(1,3); 
r21(j)= R(2,1); r22(j)= R(2,2); r23(j)= R(2,3); 
r31(j)= R(3,1); r32(j)= R(3,2); r33(j)= R(3,3); 

%end Effector Tranf Matrix
T06=[r11(j) r12(j) r13(j) px(j);
    r21(j) r22(j) r23(j) py(j);
    r31(j) r32(j) r33(j) pz(j);
    0 0 0 1];

%Inverse Kinematics Solution
%q1(j) (Hip Roll)
q1(j)= atan2(r23(j), r13(j));

if q1(j)>= pi
   
    q1(j)=0;
end

%q2(j) (Hip Yaw)
A1= [-1, -r33(j); 0, cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j)];
A2= [cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j), -1; r33(j), 0];
A= [cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j), -r33(j); r33(j), cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j)];

%other solution for Yaw q2(j) (x motion)
%A1= [-1, -r32(j); 0, -pz(j)];
%A2= [cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j), -1; cos(q1(j)).*px(j)+sin(q1(j).*py(j)), 0];
%A= [cos(q1(j)).*r13(j)+sin(q1(j)).*r23(j), -r32(j); cos(q1(j)).*px(j)+sin(q1(j)).*py(j), -pz(j)];

c2= det(A1)./det(A);
s2= det(A2)./det(A);
q2(j)= atan2(s2, c2);

if q2(j)>= pi
   
    q2(j)=0;
end

if q2(j)< 1*10^-4
    q2(j)=0;  
end

%calculate W 
T02in= [cos(q1(j)).*sin(q2(j)), sin(q1(j)).*sin(q2(j)), cos(q2(j)), 0;
    cos(q1(j)).*cos(q2(j)), cos(q2(j)).*sin(q1(j)), -sin(q2(j)), 0;
    -sin(q1(j)), cos(q1(j)), 0, 0;
    0, 0, 0, 1];

W= T02in*T06;
w11= W(1,1);
w14= W(1,4);
w34= W(3,4);
w31= W(3,1);

%q4(j)  (Knee Pitch)
x_bar= w14-l6.*w11;
y_bar= w34-l6.*w31;
c4= (x_bar.*x_bar+y_bar.*y_bar-l4.*l4-l5*l5)./(2.*l4.*l5);
s4= -sqrt(1 - c4.*c4);   %2 solutions, (knee up and knee down)
q4(j)= atan2(s4, c4);   %This finds theta2

%q3(j) (Hip Pitch)
B1= [x_bar -l5.*s4; y_bar l4+l5.*c4];
B2= [l4+l5.*c4 x_bar; l5.*s4 y_bar];
B= [l4+l5.*c4 -l5.*s4; l5.*s4 l4+l5.*c4];
c3= det(B1)./det(B);
s3= det(B2)./det(B);
q3(j)= (atan2(s3, c3));   %theta 3

%q5(j) (Ankle Pitch)
q345(j)= atan2(w31, w11);
q5(j)= (q345(j)- (q3(j)+q4(j)));    %theta 5

j=j+1;
end


%define New 5DOF DH Table
L(1)= Revolute('d',0,'a',0,'alpha',-pi/2,'offset',0);
L(2)= Revolute('d',0,'a',0,'alpha',pi/2,'offset',-pi/2);
L(3)= Revolute('d',0,'a',l4,'alpha',0,'offset',0);
L(4)= Revolute('d',0,'a',l5,'alpha',0,'offset',0);
L(5)= Revolute('d',0,'a',l6,'alpha',0,'offset',0);
bart = SerialLink(L);
%plot Bart
j=1;
for t=linspace(-T1, T1, res);

q= [q1(j) q2(j) q3(j) q4(j) q5(j)];
%bart.plot(q);
[T, ALL]= bart.fkine(q);  %Apply Forward Kinematics where ALL stores individual transformation

for i=1:5
%Converting SE3 object to a 4x4 Matrix
temp1= eye(4,3)*(tr2rt(ALL(i))*eye(3,4));
temp2= eye(4,3)*(transpose(ALL(i).transl))*eye(1,4);
temp3= circshift(temp2, [0, 3]); %Right shift 3 units
TT(:,:,i) =temp1+temp3+[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 1];   %Transforamtion Matrix for Node 1 to N, 4x4 Order
end

%cg of m4 and m5 are assumed to be at the half lengths of links
r4prime(:,1,j)= transpose((transpose(TT(:,:,3)*[-l4/2;0;0;1]))*eye(4,3));
r5prime(:,1,j)= transpose((transpose(TT(:,:,4)*[-l5/2;0;0;1]))*eye(4,3));  %contains x5prime, y5prime and z5prime 3x1 Mtx%contains x4prime, y4prime and z4prime 3x1 Mtx
r6prime(:,1,j)= transpose((transpose(TT(:,:,5)*[0;0;0;1]))*eye(4,3));
j=j+1;
end



%Now defining everything in the global coordinate/World Frame
x4= r4prime(1, :);
y4= r4prime(2, :) + yH;
z4= r4prime(3, :) + z3;

x5= r5prime(1, :);
y5= r5prime(2, :) + yH;
z5= r5prime(3, :) + z3;

x6= r6prime(1, :);
y6= r6prime(2, :) + yH;
z6= r6prime(3, :) + z3;


%Convert into Polynomials
t=linspace(-T1, T1, res);
x4_c= polyfit(t, x4, 5);
y4_c= polyfit(t, y4, 5);
z4_c= polyfit(t, z4, 5);

x5_c= polyfit(t, x5, 5);
y5_c= polyfit(t, y5, 5);
z5_c= polyfit(t, z5, 5);

x6_c= polyfit(t, x6, 5);
y6_c= polyfit(t, y6, 5);
z6_c= polyfit(t, z6, 5);

%Find the First oder erivatives
x4_dot_c= polyder(x4_c);
y4_dot_c= polyder(y4_c);
z4_dot_c= polyder(z4_c);

x5_dot_c= polyder(x5_c);
y5_dot_c= polyder(y5_c);
z5_dot_c= polyder(z5_c);

x6_dot_c= polyder(x6_c);
y6_dot_c= polyder(y6_c);
z6_dot_c= polyder(z6_c);

%Find the 2nd oder derivatives
x4_2dot_c= polyder(x4_dot_c);
y4_2dot_c= polyder(y4_dot_c);
z4_2dot_c= polyder(y4_dot_c);

x5_2dot_c= polyder(x5_dot_c);
y5_2dot_c= polyder(y5_dot_c);
z5_2dot_c= polyder(y5_dot_c);

x6_2dot_c= polyder(x6_dot_c);
y6_2dot_c= polyder(y6_dot_c);
z6_2dot_c= polyder(y6_dot_c);

%Now we need to apply the method of MMIPM and find the excitation function
z= [z4; z5; z6];
y= [y4; y5; y6];

z_dot= [polyval(z4_dot_c,t); polyval(z5_dot_c,t); polyval(z6_dot_c,t)];
y_dot= [polyval(y4_dot_c,t); polyval(y5_dot_c,t); polyval(y6_dot_c,t)];

z_2dot= [polyval(z4_2dot_c,t); polyval(z5_2dot_c,t); polyval(z6_2dot_c,t)];
y_2dot= [polyval(y4_2dot_c,t); polyval(y5_2dot_c,t); polyval(y6_2dot_c,t)];

m = [m4 m5 m6];



%% 
%Calculating the Excitation funciton for MMIPM
fMM=0;
for i=1:3
   fMM= fMM+ (m(i)./(m3_t.*yH).*(g.*z(i,:) + z(i,:).*y_2dot(i,:) - z_2dot(i,:).*y(i,:)));
end

t=linspace(-T1, T1, res);
fMM_c= polyfit(t, fMM, 5);

%Now, Solve the Diff Equ for the MMIPM Method using bvp4c
t=linspace(-T1, T1, res);
solinit = bvpinit(t, @guess);
sol = bvp4c(@(t,y) [y(2), polyval(fMM_c,t) + lamda.^2.*y(1)], @(ya, yb) [ya(1) + SW_z/2, yb(1) - SW_z/2], solinit);
z3_new= sol.y(1,:);

%convert z3_new into a polynomial
t=linspace(-T1, T1, res);
z3_new_c= polyfit(t, z3_new, 5);

%calculate error
error= max((z3_new-z3).^2)

%update the new value of z3
z3=z3_new;
z3_c= z3_new_c;
z3_dot_c= polyder(z3_c);
z3_2dot_c= polyder(z3_dot_c);

k=k-1;
end
%%
index=j;

% figure;
% t=linspace(-T1, T1, res);
% plot(t,x4);

figure;
t=linspace(-T1, T1, res);
plot(t, polyval(z3_c, t), 'r');
xlabel('t');
ylabel('z3 in [m] MMIPM');

%% Convert q into a polynomial and calculate derivatives
t=linspace(-T1, T1, res);
q1_c= polyfit(t, q1, 7);
q2_c= polyfit(t, q2, 7);
q3_c= polyfit(t, q3, 7);
q4_c= polyfit(t, q4, 7);
q5_c= polyfit(t, q5, 7);

v1_c= polyder(q1_c);
v2_c= polyder(q2_c);
v3_c= polyder(q3_c);
v4_c= polyder(q4_c);
v5_c= polyder(q5_c);

a1_c= polyder(v1_c);
a2_c= polyder(v2_c);
a3_c= polyder(v3_c);
a4_c= polyder(v4_c);
a5_c= polyder(v5_c);

%% find values of v and a now:
v1= polyval(v1_c, t);
v2= polyval(v2_c, t);
v3= polyval(v3_c, t);
v4= polyval(v4_c, t);
v5= polyval(v5_c, t);

a1= polyval(a1_c, t);
a2= polyval(a2_c, t);
a3= polyval(a3_c, t);
a4= polyval(a4_c, t);
a5= polyval(a5_c, t);


%% For complete robot walk. Left and right both
t=linspace(0, 4*T1, 2*res);
for j=1:index
    
q1(index-2+j)= q1(index-1);
q2(index-2+j)= q2(index-1);
q3(index-2+j)= q3(index-1);
q4(index-2+j)= q4(index-1);
q5(index-2+j)= q5(index-1);

v1(index-2+j)= v1(index-1);
v2(index-2+j)= v2(index-1);
v3(index-2+j)= v3(index-1);
v4(index-2+j)= v4(index-1);
v5(index-2+j)= v5(index-1);

a1(index-2+j)= a1(index-1);
a2(index-2+j)= a2(index-1);
a3(index-2+j)= a3(index-1);
a4(index-2+j)= a4(index-1);
a5(index-2+j)= a5(index-1);

end

t=linspace(0, 4*T1, 2*res);
for j=1:index-1
    q3R(j)= q3(2*(index-1));
    q4R(j)= q4(2*(index-1));
    q5R(j)= q5(2*(index-1));
end

for j=1:index-1
   q3R(index-1+j)= q3(j);
   q4R(index-1+j)= q4(j);
   q5R(index-1+j)= q5(j);
    
end


if ~exist('actuatorType','var')
    actuatorType=1
    
end

%% Foot Dimension
foot_x=0.09;  %9 cm
foot_y=0.05;
foot_z=0.01;

torso_x=0.1;
torso_y=0.1;
torso_z=0.15;


%% Inertia Matrices for each link
I1= zeros(3,3);
I2= zeros(3,3);
I3= 10^-9*[358897 -0.38 0.26;
    -0.38 133672 -162883;
    0.26 -162883 273096];
I4= 10^-9*[271290 0.40 0.16;
    0.40 49662 89240;
    0.16 89240 248535];
I5= 10^-9*[19353 74 2.5;
    74 29106 84;
    2.5 84 22353];

%cg poistions for each link 3x1 matrix
r1=[0; 0; 0];
r2=[0; 0; 0];
r3=[-l4/2; 0; 0];
r4=[-l5/2; 0; 0];
r5=[0; 0; 0];

%gains and offsets for Servo Control
q2L_offset=0;
q3L_offset=94;
q4L_offset=-29;
q5L_offset=-10;

q2R_offset=0;
q3R_offset=70;
q4R_offset=-8;
q5R_offset=1;

%gains
q4R_gain=-1;


function g = guess(t)
g = [0
     0];
end

function output3= RX(cc) 
output3= [1 0 0;
    0 cos(cc) -sin(cc);
    0 sin(cc) cos(cc)];
end

function output2= RY(bb)
output2= [cos(bb) 0 sin(bb);
    0 1 0;
    -sin(bb) 0 cos(bb)];
end

function output1= RZ(aa)
output1= [cos(aa) -sin(aa) 0;
    sin(aa) cos(aa) 0;
    0 0 1];
end




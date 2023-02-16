clear;
clc;
close all;
hold on;

i = 2;
view(3);
set(gca,'DataAspectRatio',[1 1 1]);
Ni = 50;
len1 = 150e-3;
global R10 R20 ez curvature1 curvature2 alpha
alpha = 180; %Base rotation of tube 1 and tube 2
R10 = eye(3);%Frame rotation matrix of tube 1
R20 = eye(3)*rotz(alpha);%Frame rotation matrix of tube 2
ez = [0 0 1]';
global v k1xy k2xy k1z k2z K1 K2
v = 0.3;

k1xy = 5.07e-2; k2xy = 5.07e-2;
k1z = k1xy/(1+v);k2z = k2xy/(1+v);

K1 = diag([k1xy, k1xy, k1z]);
K2 = diag([k2xy, k2xy, k2z]);

len2 = 150e-3;
curvature1 = 1/(150e-3);
curvature2 = 1/(150e-3);

global u1_pre_const u2_pre_const
u1_pre = @(s) [0 1 0]*curvature1;
u2_pre = @(s) [0 1 0]*curvature2;
u1_pre_const = u1_pre(1);
u2_pre_const = u2_pre(1);

uhat_1 = zeros(3*Ni,1);% Inner tube center line precurvature
phat_1 = zeros(3*Ni,1);% Inner tube center line position
for m = 1:Ni
    uhat_1(3*m-2:3*m) = u1_pre(1);
end
[phat_1,~] = findShape2(uhat_1,len1,Ni,R10,"Tube 1");



uhat_2 = zeros(3*Ni,1);% Outer tube center line precurvature
phat_2 = zeros(3*Ni,1);% Outer tube center line position
for m = 1:Ni
    uhat_2(3*m-2:3*m) = u2_pre(1);
end
[phat_2,~] = findShape2(uhat_2,len2,Ni,R20,"Tube 2");

ini = [1 90 reshape(R20,[1,9]) 0 0 0];%[u2z_ini,a2_ini,R_ini,p_ini,u1z_ini]
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
best_ini = fsolve(@residual_fun,ini);
[s y] = ode45(@twoTube,[0 len2],best_ini);

u2z_array = reshape(y(:,1)',[],1);
p_array_2 = y(:,12:14);
plot3(p_array_2(:,1),p_array_2(:,2),p_array_2(:,3),"DisplayName","Final Solution From p");

%Recover the whole u1,u2 vector
u2 = [];
u1 = [];
for i = 1:size(s)
    u1_hat = u1_pre(1);
    u2_hat = u2_pre(1);
    a1 = 0;
    a2 = y(i,2);
    Rz_a1 = rotz(a1);
    Rz_a2 = rotz(a2);
    M = Rz_a1*K1*u1_hat' + Rz_a2*K2*u2_hat';
    u2xy_s = (inv(K1+K2)*Rz_a2'*M);
    u1xy_s = (inv(K1+K2)*Rz_a1'*M);
    u2 = [u2;u2xy_s(1:2);y(i,1)];
    u1 = [u1;u1xy_s(1:2);(-1/k1z)*(k2z*y(i,1))]
end

g2 = zeros(4,4,size(s,1));
for i = 1:size(s)
    R = reshape(y(i,3:11),3,3);
    p = p_array_2(i,:);
    g2(:,:,i) = [R p';0 0 0 1];
end

% Use findShape2 to validate if it is correctly written
findShape2(u2,len1,size(y,1),R20,"Final Solution From u");
legend();
legend show;

save("ini.mat","u2","u1");
%save("zzy.mat","p_array_2");
function dyds = twoTube(s,y)
    %y = [u2z;a2;R;p]    
    global v k1xy k2xy k1z k2z K1 K2
    
    u2z = y(1);
    a2 = y(2); %Rotation from tube one
    a1 = 0;    %Rotation from tube one, 0 by definition
    R = reshape(y(3:11),[3,3]);
    p = y(12:14);

    global ez u1_pre_const u2_pre_const
    
    Rz_a1 = rotz(a1);
    Rz_a2 = rotz(a2);
    
    M = Rz_a1*K1*u1_pre_const' + Rz_a2*K2*u2_pre_const';
    u2xy_s = (Rz_a2'*inv(K1+K2)*M);
    u1z = (-1/k1z)*(k2z*u2z);

    du2z_ds = (k2xy/k2z)*(u2xy_s(1)*u2_pre_const(2)-u2xy_s(2)*u2_pre_const(1));
    
    u2 = [u2xy_s(1:2);u2z];

    %u2_frame0 = Rz_a2*u2;
    a2_dot = u2z-u1z;

    dR_ds = R*hat(u2);
    dp_ds = R*ez;
    dyds = [du2z_ds;a2_dot;reshape(dR_ds,[9,1]);dp_ds];
end

function res = residual_fun(ini)
    %ini is a initial guess which will be used for ode solving
    global alpha
    [~,sol] = ode45(@twoTube,[0,150e-3],ini);
    
    u2z = sol(:,1);
    a2 = sol(:,2);
    res = zeros(2,1);

    res(1) = (u2z(end) - 0);%namely uiz(L) = 0
    res(2) = (a2(1) - alpha); %namely a2(0) = theta2(0) - theta1(0)
    %res(3) = (u1z(end) - 0);
    %disp(res);
    %disp(res);
end
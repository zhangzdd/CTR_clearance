clear;
clc;
close all;
hold on;

tubeNumber = 2;
view(3);
set(gca,'DataAspectRatio',[1 1 1]);
Ni = 50;
len1 = 150e-3;
len2 = 150e-3;
alpha = 180;
curvature1 = 1/(150e-3);
curvature2 = 1/(150e-3);

u1_pre = @(s) [0 1 0]*curvature1;
u2_pre = @(s) [0 1 0]*curvature2;

tube1 = tube(150e-3,Ni,u1_pre,[20.07e-2,0.3]);
tube2 = tube(150e-3,Ni,u1_pre,[5.07e-2,0.3]);
tset = tubeset({tube1,tube2},{[0,0],[180,0]});
R10 = eye(3);%Frame rotation matrix of tube 1
R20 = eye(3)*rotz(alpha);%Frame rotation matrix of tube 2
ez = [0 0 1]';

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

ini = [1;90;reshape(R20,[9,1]);0;0;0];%[u2z_ini,a2_ini,R_ini,p_ini,u1z_ini]
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
f = @(y) boundary_condition(y,tset); % function of 'real' variable y,tset embedded as constant
best_ini = fsolve(f,ini);

[u2out,a2out,pout] = integrate_tube(best_ini,tset);

u2z_array = reshape(u2out,[],1);
p_array_2 = pout;
plot3(p_array_2(1,:),p_array_2(2,:),p_array_2(3,:),"DisplayName","Final Solution From p");

%Recover the whole u1,u2 vector
u2 = [];
u1 = [];
for i = 1:Ni
    u1_hat = u1_pre(1);
    u2_hat = u2_pre(1);
    a1 = 0;
    a2 = a2out(i);
    Rz_a1 = rotz(a1);
    Rz_a2 = rotz(a2);
    M = Rz_a1*tube1.K*u1_hat' + Rz_a2*tube2.K*u2_hat';
    u2xy_s = (inv(tube1.K+tube2.K)*Rz_a2'*M);
    u1xy_s = (inv(tube1.K+tube2.K)*Rz_a1'*M);
    u2 = [u2;u2xy_s(1:2);u2out(3,i)];
    u1 = [u1;u1xy_s(1:2);(-1/tube1.kz)*(tube2.kz*u2out(3,i))]
end

% g2 = zeros(4,4,size(s,1));
% for i = 1:size(s)
%     R = reshape(y(i,3:11),3,3);
%     p = p_array_2(i,:);
%     g2(:,:,i) = [R p';0 0 0 1];
% end

% Use findShape2 to validate if it is correctly written
findShape2(u2,len1,Ni,R20,"Final Solution From u");
legend();
legend show;
save("ini.mat","u2","u1");
%warning('off','last');
close all;
clear;
clc;
hold on;
set(gca,'DataAspectRatio',[1 1 1]);
view(3);

legend show;
i = 2; %number of tubes 

%Load zero clearance solution (initial guess) and discretization points
ini = load("ini.mat");
alpha1 = ini.u1;
alpha2 = ini.u2;

Ni = size(alpha1,1)/3; % Discretize coefficeint, based on zero clearance solution 
N = i*Ni;
ez = [0 0 1];
len1 = 500e-3;

curvature1 = 1/(200e-3);
curvature2 = 1/(250e-3);

len2 = len1;
delta_s = len1/(Ni-1);
delta = 10;


R10 = eye(3);%Base rotation of tube 1
R20 = eye(3)*expm(hat([0 0 1]*pi));%Base rotation of tube 2


uhat_1 = zeros(3*Ni,1);% Inner tube center line precurvature
phat_1 = zeros(3*Ni,1);% Inner tube center line position

for m = 1:Ni
    uhat_1(3*m-2:3*m) = [0 1 0]*curvature1;%Local frame curvature
end
[phat_1,R1_hat,g1_hat] = findShape2(uhat_1,len1,Ni,R10,"Tube 1");

uhat_2 = zeros(3*Ni,1);% Outer tube center line precurvature
phat_2 = zeros(3*Ni,1);% Outer tube center line position

% uhat_2(1:3*floor(Ni/2)) = repmat([0 0 0],1,floor(Ni/2));
% uhat_2(3*floor(Ni/2)+1:3*floor(Ni/2)+3) = [0 1/(delta_s) 0];
% uhat_2(3*floor(Ni/2)+4:3*Ni) = repmat([0 0 0],1,floor(Ni/2));
% phat_2 = findShape2(uhat_2,len2,Ni,R20,"Tube 2");

% This part is for the tunnel project

for m = 1:Ni
    uhat_2(3*m-2:3*m) = [0 1 0]*curvature2;
end
[phat_2,R2_hat,g2_hat] = findShape2(uhat_2,len2,Ni,R20,"Tube 2");

uhat = [uhat_1;uhat_2];


v1 = 0.3; %Poisson Ratio
k1xy = 5.07e-2;k1z = k1xy/(1+v1);% Bending and torsional stiffness of inner tube
k = [k1xy k1xy k1z];Kd = k;
for m = 1:Ni-1
    Kd = [Kd k];
end

v2 = 0.3;
k2xy = 20.07e-2;k2z = k2xy/(1+v2);% Bending and torsional stiffness of inner tube
k = [k2xy k2xy k2z];
for m = 1:Ni
    Kd = [Kd k];
end
Kd = Kd*delta_s;
%K = diag(Kd);
K = diag(Kd);

c = 200e-3;% 3mm clearance
max_clearance_step = 10;



% Calculate S
S_zzy = zeros(Ni,3,3*Ni*i);
for m = 1:Ni
    %S(m,n,:,:) = zeros(3,3*Ni*i); Shape indications
    S_zzy(m,:,3*m-2:3*m) = eye(3);
    S_zzy(m,:,3*Ni+3*m-2:3*Ni+3*m) = -eye(3);
end

%Outer loop for optimizing u,p


%I tried to add some deviation but that worked awfully
ustar_1_deviation = 0*(rand(3*Ni,1)-1/2);
ustar_2_deviation = 0*(rand(3*Ni,1)-1/2);


% Start the initial guess with zero clearance solutions
ustar_1 = reshape(alpha1,3*Ni,1)+ustar_1_deviation;
ustar_2 = reshape(alpha2,3*Ni,1)+ustar_2_deviation;

[pstar_1,Rstar_1,g1_star] = findShape2(ustar_1,len1,Ni,R10,"Tube 1 initial guess");
[pstar_2,Rstar_2,g2_star] = findShape2(ustar_2,len2,Ni,R20,"Tube 2 initial guess");
[p_alpha,~,g_alpha] = findShape2(reshape(alpha2,3*Ni,1),len1,Ni,R20,"Zero Clearance Solution");


%disp(norm(pstar_1 - pstar_2));

ustar = [ustar_1;ustar_2];
pstar = [pstar_1;pstar_2];

count = 0;

for clearance_step = 1:max_clearance_step
% Loop for updating guess of u

%Update ustar1/2 from updated ustar
ustar_1 = ustar(1:3*Ni);
ustar_2 = ustar(3*Ni+1:i*Ni*3);
pstar_1 = pstar(1:3*Ni);
pstar_2 = pstar(3*Ni+1:i*Ni*3);

g = K*(ustar - uhat);

% Calculate Jp
Jp1 = zeros(3*Ni,3*Ni);
Jp2 = zeros(3*Ni,3*Ni);
% m,n for j,k

for m = 1:Ni
    for n = 1:Ni
        if m<=n
            Jp1(3*m-2:3*m,3*n-2:3*n) = zeros(3,3);
            continue;
        end
        %Rstar_k = expm(hat(ustar_1(3*n-2:3*n)));
        Rstar_k = squeeze(Rstar_1(n,:,:));
        Jp1(3*m-2:3*m,3*n-2:3*n)=hat(pstar_1(3*n-2:3*n)-pstar_1(3*m-2:3*m))*Rstar_k*delta_s;
    end
end

for m = 1:Ni
    for n = 1:Ni
        if m<=n
            Jp2(3*m-2:3*m,3*n-2:3*n) = zeros(3,3);
            continue;
        end
        %Rstar_k = expm(hat(ustar_2(3*n-2:3*n)));
        Rstar_k = squeeze(Rstar_2(n,:,:));
        Jp2(3*m-2:3*m,3*n-2:3*n)=hat(pstar_2(3*n-2:3*n)-pstar_2(3*m-2:3*m))*Rstar_k*delta_s;
    end
end
Jp = blkdiag(Jp1,Jp2);

%Calculate P and pdot
P1 = zeros(Ni,3,3);  
P2 = zeros(Ni,3,3);
pdot_1_array = {};
pdot_2_array = {};
for j = 1:Ni
    pdot_1 = squeeze(Rstar_1(j,:,:))*ez';
    pdot_1_array{j} = pdot_1;
    P1(j,:,:) = eye(3)-pdot_1*pdot_1';
end 

for j = 1:Ni
    pdot_2 = squeeze(Rstar_2(j,:,:))*ez';
    pdot_2_array{j} = pdot_2;
    P2(j,:,:) = eye(3)-pdot_2*pdot_2';
end
P_zzy = zeros(i,Ni,3,3);
P_zzy(1,:,:,:) = P1;
P_zzy(2,:,:,:) = P2;


%Calculate X, with respect to tube 1
X_zzy = zeros(Ni,3,3*Ni*i);
for n = 1:Ni
    X_zzy(n,:,:) = squeeze(P_zzy(2,n,:,:))*squeeze(S_zzy(n,:,:))*Jp;
end


%Calculate q
q = zeros(i-1,Ni);
for m = 1:i-1
    for n = 1:Ni
        q(m,n) = 1/2*(c/max_clearance_step)^2;
    end
end

% Formulate h and q from lambda and optimize upon lambda
%lambda = ones((i-1)*Ni,1);
previous = 0;
previous_lambda = zeros((i-1)*Ni,1);
maxIter = 200;
stepSize = 0.5;
lambda_zzy = 0.1/Ni*ones(Ni,1);
dlambda_zzy = zeros(Ni,1);
%lambda optimizing loop
for lp = 1:maxIter
    
    %disp(lp);

    relativeSize = stepSize*(norm(lambda_zzy)+1e-6);
    dlambda_zzy = relativeSize * dlambda_zzy/(norm(dlambda_zzy)+1e-9);
    lambda_zzy = lambda_zzy + dlambda_zzy;
    lambda_zzy(lambda_zzy<0) = 0;
    
    Q_lambda = K;
    for m = 1:i-1
        for n = 1:Ni
         k = (m-1)*Ni+n;
         Q_lambda = Q_lambda + lambda_zzy(k)*squeeze(X_zzy(n,:,:))'*squeeze(X_zzy(n,:,:));
        end
    end

    h_lambda = g;
    delta_u_zzy =  -Q_lambda\h_lambda;

    for m = 1:i-1
        for n = 1:Ni
        k = (m-1)*Ni+n;
        % Compute the angle between two vectors
%         pdot_1 = pdot_1_array{n};
%         pdot_2 = pdot_2_array{n};
%         angle = acos(dot(pdot_1,pdot_2)/(norm(pdot_1)*norm(pdot_2)));
%         
        G_u = 1/2*delta_u_zzy'*squeeze(X_zzy(n,:,:))'*squeeze(X_zzy(n,:,:))*delta_u_zzy-q(m,n);
        dlambda_zzy(k) = G_u;                        
        end
    end
    
    
    q_col = reshape(q,(i-1)*Ni,1);
    J_lambda_zzy = -1/2*g'*inv(Q_lambda)*g-q_col'*lambda_zzy;
    %disp(J_lambda_zzy-previous)
    disp(J_lambda_zzy)
% 
%     if abs(J_lambda_zzy-previous)<9e-6 
% %         if J_lambda>previous
% %             break;
% %         end
%         break;
%         
%     end
    
    
    if J_lambda_zzy < previous
        if lp>1
        stepSize = stepSize*0.5;
        end
    end
    previous = J_lambda_zzy;
    figure(5);
    plot(lambda_zzy);
    end
    
    
    % lambda already optimized, use current delta_u for update
    
    
    %obj = 1/2*delta_u'*K*delta_u+g'*delta_u;
    count = count+1;
    ustar = ustar+delta_u_zzy;
    pstar = pstar+Jp*delta_u_zzy;


    [~,R1_new,g1_new] = findShape2(ustar(1:3*Ni),len1,Ni,R10,['Tube 1, iter ',num2str(count)]);
    [~,R2_new,g2_new] = findShape2(ustar(3*Ni+1:3*i*Ni),len1,Ni,R20,['Tube 2, iter ',num2str(count)]);
    pstar_1_new = pstar(1:3*Ni);
    pstar_2_new = pstar(3*Ni+1:3*i*Ni);
    
    Rstar_2 = R2_new;
    Rstar_1 = R1_new;
    
    %for each point on tube 1 center line, search for the nearest point in
    %tube 2 center line
    tube1_pointset = reshape(pstar_1_new,3,Ni)';
    tube2_pointset = reshape(pstar_2_new,3,Ni)';
    %create a new mapping to map segment in tube 1 to segment in tube 2
    d = dictionary;
    for m = 1:Ni
        Idx = knnsearch(tube2_pointset,tube1_pointset(m,:));
        d(m) = Idx;
    end
    %Update S based on new mapping
    S_zzy = zeros(Ni,3,3*Ni*i);
    for m = 1:Ni
        %S(m,n,:,:) = zeros(3,3*Ni*i); Shape indications
        S_zzy(m,:,3*m-2:3*m) = eye(3);
        S_zzy(m,:,3*Ni+3*d(m)-2:3*Ni+3*d(m)) = -eye(3);
    end
    %Make plot
    h = figure(2);

    rin1 = 0.95e-3;
    rout1 = 1.15e-3;

    rin2 = rout1+(c/max_clearance_step*clearance_step);
    rout2 = rin2 + 2e-3;
%     
    hold off;
    plot3DTubes(g1_new, rin1, rout1, 1,[1,0,0]);
    hold on;
    plot3DTubes(g2_new, rin2, rout2, 0.15);
    %plot3DTubes(g1_hat,rin1, rout1, 1,[1,0,0]);
    title(num2str(clearance_step));
    view(0,0);
    axis equal;
    %save("p_zzy.mat","pstar_2_new","pstar_1_new");
end
gstar_1_new_lc = g1_new;
gstar_2_new_lc = g2_new;
save("pstar_large_clearance.mat","gstar_1_new_lc","gstar_2_new_lc");

rin1 = 0.95e-3;
rout1 = 1.15e-3;

rin2 = rout1+c;
rout2 = rin2 + 2e-3;
figure(3);subplot(1,3,1);
plot3DTubes(g1_hat, rin1, rout1, 0.1,[1,0,0]);
axis equal;
view(0,0);

subplot(1,3,2);
plot3DTubes(g2_hat, rin2, rout2, 0.15);
axis equal;
view(0,0);

subplot(1,3,3);
plot3DTubes(g1_new, rin1, rout1, 1,[1,0,0]);
hold on;
plot3DTubes(g2_new, rin2, rout2, 0.15);
axis equal;
view(0,0);

%figure(2);axis([-1 1 -1 1 -1 1]);xlabel('X');ylabel('Y');zlabel('Z');view(0,90);pause;view(0,0);pause;view(90,0);
disp("Finished")
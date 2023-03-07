clear;
clc;
close all;
clear;
clc;
hold on;
set(gca,'DataAspectRatio',[1 1 1]);
view(3);
%parpool(4);
legend show;

r1out = 1.15e-3;
r2in = 2.15e-3;
r1thickness = 0.2e-3;
r2thickness = 2e-3;

circle_sample = 4;
theta = linspace(0,2*pi,circle_sample);
c1 = {};
for t = theta
    p = [r1out*cos(t) r1out*sin(t) 0]';
    c1{end+1} = p;
end

i = 2; %number of tubes 

%Load zero clearance solution (initial guess) and discretization points
ini = load("ini.mat");
alpha1 = ini.u1;
alpha2 = ini.u2;

Ni = size(alpha1,1)/3; % Discretize coefficeint, based on zero clearance solution 
N = i*Ni;
ez = [0 0 1];
len1 = 150e-3;

curvature1 = 1/(50e-3);
curvature2 = 1/(150e-3);

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

%This part is for the tunnel project

for m = 1:Ni
    uhat_2(3*m-2:3*m) = [0 1 0]*curvature2;
end
[phat_2,R2_hat,g2_hat] = findShape2(uhat_2,len2,Ni,R20,"Tube 2");

uhat = [uhat_1;uhat_2];


v1 = 0.3; %Poisson Ratio
k1xy = 20.07e-2;k1z = k1xy/(1+v1);% Bending and torsional stiffness of inner tube
k = [k1xy k1xy k1z];Kd = k;
for m = 1:Ni-1
    Kd = [Kd k];
end

v2 = 0.3;
k2xy = 5.07e-2;k2z = k2xy/(1+v2);% Bending and torsional stiffness of inner tube
k = [k2xy k2xy k2z];
for m = 1:Ni
    Kd = [Kd k];
end
Kd = Kd*delta_s;
%K = diag(Kd);
K = diag(Kd);

%c = 1e-3;% 3mm clearance
max_clearance_step = 1;



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
% Loop for updating guess

while 1

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

%Calculate P 
P1 = zeros(Ni,3,3);  
P2 = zeros(Ni,3,3);
for j = 1:Ni
    pdot_1 = squeeze(Rstar_1(j,:,:))*ez';
    P1(j,:,:) = eye(3)-pdot_1*pdot_1';
end 

for j = 1:Ni
    pdot_2 = squeeze(Rstar_2(j,:,:))*ez';
    P2(j,:,:) = eye(3)-pdot_2*pdot_2';
end
P_zzy = zeros(i,Ni,3,3);
P_zzy(1,:,:,:) = P1;
P_zzy(2,:,:,:) = P2;

% Calculate S
S_zzy = zeros(1,Ni,3,3*Ni*i);
S1 = repmat(eye(3),1,Ni);
S2 = repmat(-eye(3),1,Ni);
for m = 1:i-1
    for n = 1:Ni
        %S(m,n,:,:) = zeros(3,3*Ni*i); each element should be a 3 by 3*N matrix
        %S_zzy(m,n,:,3*Ni*(m-1)+3*n-2:3*Ni*(m-1)+3*n) = squeeze(P2(n,:,:));
        S_zzy(m,n,:,3*Ni*(m-1)+3*n-2:3*Ni*(m-1)+3*n) = eye(3);
        S_zzy(m,n,:,3*Ni*m+3*n-2:3*Ni*m+3*n) = -eye(3);
    end
end




% Formulate h and q from lambda and optimize w.r.t lambda
%lambda = ones((i-1)*Ni,1);
previous = 0;
previous_lambda = zeros((i-1)*Ni*circle_sample,1);
maxIter = 200;

lambda_zzy = 0.1/Ni*ones(Ni*circle_sample,1);
stepSize = 0.5;

dlambda_zzy = zeros(Ni*circle_sample,1);
%lambda optimizing loop
for lp = 1:maxIter
    %disp(lp);
    relativeSize = stepSize*(norm(lambda_zzy)+1e-6);
    dlambda_zzy = relativeSize * dlambda_zzy/(norm(dlambda_zzy)+1e-9);
    lambda_zzy = lambda_zzy + dlambda_zzy;
    lambda_zzy(lambda_zzy<0) = 0;
    
    Q_lambda = K;
    % Here m is the label for ci
    for m = 1:circle_sample
        for n = 1:Ni
         k = (m-1)*Ni+n;
         Sj = squeeze(S_zzy(1,n,:,:));
         P2j = squeeze(P2(n,:,:));
         Q_lambda = Q_lambda + lambda_zzy(k)*2*Jp'*(Sj'*P2j*Sj)*Jp;
        end
    end

    h_lambda = g;        
    for n = 1:Ni
         Sj = squeeze(S_zzy(1,n,:,:));
         R1j = squeeze(Rstar_1(n,:,:));
         P2j = squeeze(P2(n,:,:));
        for m = 1:circle_sample
         ci = c1{m};
         k = (m-1)*Ni+n;
         h_lambda = h_lambda + lambda_zzy(k)*2*((ci'*R1j'*P2j*Sj*Jp)'+(pstar'*(Sj'*P2j*Sj)*Jp)');
        end
    end

    invQ = inv(Q_lambda);
    delta_u_zzy = -invQ*h_lambda;

    for n = 1:Ni
        Sj = squeeze(S_zzy(1,n,:,:));
        R1j = squeeze(Rstar_1(n,:,:));
        P2j = squeeze(P2(n,:,:));
        dh_dlambdak = 2*((ci'*R1j'*P2j*Sj*Jp)'+(pstar'*(Sj'*P2j*Sj)*Jp)');
        dQ_dlambdak = (2*Jp'*Sj'*P2j*Sj*Jp);
        dinvQ_dlambdak = -invQ*dQ_dlambdak*invQ;
        dH1_dlambdak = 2*h_lambda'*invQ*dh_dlambdak + h_lambda'*dinvQ_dlambdak*h_lambda;
        for m = 1:circle_sample   
            ci = c1{m};
            k = (m-1)*Ni+n;
            %H = H1+H2          
            dH2_dlambdak = ci'*R1j'*P2j*R1j*ci + 2*ci'*R1j'*P2j*Sj*pstar + pstar'*Sj'*P2j*Sj*pstar-r2in^2;
            dlambda_zzy(k) = -dH1_dlambdak+dH2_dlambdak;                        
        end
    end

    %disp(mean(dlambda_zzy));
    J_lambda_zzy = -1/2*h_lambda'*invQ*h_lambda+compute_coef_sum_V(c1,Rstar_1,P2,Jp,S_zzy,pstar,Ni,circle_sample,lambda_zzy)-sum(lambda_zzy.*r2in^2);
    %disp(J_lambda_zzy)

    
%     if mod(lp,50)==0 && lp > 2
%         stepSize = stepSize*0.6;
%     end
    
    if norm(J_lambda_zzy - previous)<1e-7
        break;
    end
    
    if J_lambda_zzy < previous
        if lp>1
           stepSize = stepSize*0.5;
        %break;
        end
    end
    previous = J_lambda_zzy;

   end
    
    
    % lambda already optimized, use current delta_u for update
    
    
    %flag = check_penetration(c1,Rstar_1,P2,Jp,S_zzy,pstar,Ni,circle_sample,lambda_zzy,delta_u_zzy,r2in);
    %disp(flag);
    count = count+1;
    ustar = ustar+delta_u_zzy;

    [pstar_1_new,R1_new,g1_new] = findShape2(ustar(1:3*Ni),len1,Ni,R10,['Tube 1, iter ',num2str(clearance_step)]);
    [pstar_2_new,R2_new,g2_new] = findShape2(ustar(3*Ni+1:3*i*Ni),len1,Ni,R20,['Tube 2, iter ',num2str(clearance_step)]);
    h = figure(2);
    set(h,'Position', [100,100, 500, 850]);
    view([0,0]);

    hold off
    plot3DTubes(g1_new, r1out-r1thickness, r1out, 1,[1,0,0]);
    hold on
    plot3DTubes(g2_new, r2in, r2in + r2thickness, 0.15);
    plot3DTubes(g1_hat,r1out-r1thickness, r1out, 1,[1,0,0]);
    axis equal;
    view(0,0);
    
    obj = 1/2*(ustar-uhat)'*K*(ustar-uhat);
    disp(['Elastic energy ',num2str(obj)]);
    disp(norm(delta_u_zzy))
    gstar_1_new_ellip = g1_new;
    gstar_2_new_ellip = g2_new;
    save("pstar_ellip.mat","gstar_1_new_ellip","gstar_2_new_ellip");

    %plotSection(g2_new,3e-3,[0 1 0]);
    %Compare updating with Jp and updating with findShape
    
    Rstar_2 = R2_new;
    Rstar_1 = R1_new;
    pstar = [pstar_1_new; pstar_2_new];
    %pstar = pstar_Jp;
    %lambda_from_last_guess = lambda_zzy;
    %stepSize_from_last_guess = stepSize;
end
end





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
alpha = ini.u2;

Ni = size(alpha,1)/3; % Discretize coefficeint, based on zero clearance solution 
N = i*Ni;
ez = [0 0 1];
len1 = 150e-3;

curvature1 = 1/(150e-3);
curvature2 = 1/(150e-3);

len2 = len1;
delta_s = len1/Ni;
delta = 1e-6;


R10 = eye(3);%Base rotation of tube 1
R20 = eye(3)*expm(hat([0 0 1]*pi));%Base rotation of tube 2


uhat_1 = zeros(3*Ni,1);% Inner tube center line precurvature
phat_1 = zeros(3*Ni,1);% Inner tube center line position

for m = 1:Ni
    uhat_1(3*m-2:3*m) = [1 0 0]*curvature1;%Local frame curvature
end
[phat_1,R1_hat,~] = findShape2(uhat_1,len1,Ni,R10,"Tube 1");

uhat_2 = zeros(3*Ni,1);% Outer tube center line precurvature
phat_2 = zeros(3*Ni,1);% Outer tube center line position

% uhat_2(1:3*Ni/2) = repmat([0 0 0],1,Ni/2);
% uhat_2(3*Ni/2+1:3*Ni) = repmat([-pi/2 0 0],1,Ni/2);
% phat_2 = findShape2(uhat_2,len2,Ni);
% This part is for the tunnel project

for m = 1:Ni
    uhat_2(3*m-2:3*m) = [1 0 0]*curvature2;
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
k2xy = 5.07e-2;k2z = k2xy/(1+v2);% Bending and torsional stiffness of inner tube
k = [k2xy k2xy k2z];
for m = 1:Ni
    Kd = [Kd k];
end

K = diag(Kd);

c = 3e-3;% 3mm clearance



% Calculate S
S = zeros(1,Ni,3,3*Ni*i);
S1 = repmat(eye(3),1,Ni);
S2 = repmat(-eye(3),1,Ni);
for m = 1:i-1
    for n = 1:Ni
        %S(m,n,:,:) = zeros(3,3*Ni*i); Shape indications
        S(m,n,:,3*Ni*(m-1)+3*n-2:3*Ni*(m-1)+3*n) = eye(3);
        S(m,n,:,3*Ni*m+3*n-2:3*Ni*m+3*n) = -eye(3);
    end
end

%Outer loop for optimizing u,p


%I tried to add some deviation but that worked awfully
ustar_1_deviation = 0*(rand(3*Ni,1)-1/2);
ustar_2_deviation = 0*(rand(3*Ni,1)-1/2);


% Start the initial guess with zero clearance solutions
ustar_1 = alpha+ustar_1_deviation;
ustar_2 = alpha+ustar_2_deviation;

[pstar_1,Rstar_1,g1_star] = findShape2(ustar_1,len1,Ni,R20,"Tube 1 initial guess");
[pstar_2,Rstar_2,g2_star] = findShape2(ustar_2,len2,Ni,R20,"Tube 2 initial guess");
[p_alpha,~,g_alpha] = findShape2(alpha,len1,Ni,R20,"Zero Clearance Solution");
%plotSection(g_alpha,3e-3,[0 1 0]);
%plotSection(g2_star,3e-3,[0 1 0]);


%disp(norm(pstar_1 - pstar_2));

ustar = [ustar_1;ustar_2];
pstar = [pstar_1;pstar_2];

count = 0;
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
        Jp1(3*m-2:3*m,3*n-2:3*n)=hat(pstar_1(3*n-2:3*n)-pstar_1(3*m-2:3*m))*Rstar_k;
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
        Jp2(3*m-2:3*m,3*n-2:3*n)=hat(pstar_2(3*n-2:3*n)-pstar_2(3*m-2:3*m))*Rstar_k;
    end
end
Jp = blkdiag(Jp1,Jp2);

%Calculate P 
P1 = zeros(Ni,3,3);  
P2 = zeros(Ni,3,3);
P = zeros(i,Ni,3,3);
for j = 1:Ni
    pdot_1 = squeeze(Rstar_1(j,:,:))*ez';
    P1(j,:,:) = eye(3)-pdot_1*pdot_1';
end 

for j = 1:Ni
    pdot_2 = squeeze(Rstar_2(j,:,:))*ez';
    P2(j,:,:) = eye(3)-pdot_2*pdot_2';
end
P(1,:,:,:) = P1;
P(2,:,:,:) = P2;


%Calculate X
X = zeros(i,Ni,3,3*Ni*i);
for m = 1:i-1
    for n = 1:Ni
        X(m,n,:,:) = squeeze(P(m,n,:,:))*squeeze(S(m,n,:,:))*Jp;
        %X(isnan(S)) = 0;
    end
end
       



%Calculate q
q = zeros(i-1,Ni);
for m = 1:i-1
    for n = 1:Ni
        q(m,n) = 1/2*c^2;
    end
end

% Formulate h and q from lambda and optimize upon lambda
lambda = ones((i-1)*Ni,1);
previous = 0;
previous_lambda = zeros((i-1)*Ni,1);
%lambda optimizing loop
while 1
    %Formulate Q_lambda
    Q_lambda = K;
    for m = 1:i-1
        for n = 1:Ni
         k = (m-1)*Ni+n;
         Q_lambda = Q_lambda + lambda(k)*squeeze(X(m,n,:,:))'*squeeze(X(m,n,:,:));
        end
    end

    h_lambda = g;
    delta_u = Q_lambda\h_lambda;
    for m = 1:i-1
        for n = 1:Ni
        k = (m-1)*Ni+n;
        G_u = 1/2*g'*inv(Q_lambda)*squeeze(X(m,n,:,:))'*squeeze(X(m,n,:,:))*inv(Q_lambda)*g-q(m,n);
                
                
        lambda(k) = lambda(k) + delta*G_u;
        end
   
    end
    
    
    q_col = reshape(q,(i-1)*Ni,1);
    J_lambda = 1/2*g'*inv(Q_lambda)*g-q_col'*lambda;

    disp(J_lambda);
    lambda(lambda<0) = 0;

    disp(J_lambda>previous);
    if J_lambda-previous<1e-9 
%         if J_lambda>previous
%             break;
%         end
        break;
        
    end
    if J_lambda<previous
        delta = delta*0.97;
    end
    previous = J_lambda;

    end
    
    
    % lambda already optimized, use current delta_u for update
    

    obj = 1/2*delta_u'*K*delta_u+g'*delta_u;


    count = count+1;
    ustar = ustar+delta_u;
    pstar_Jp = Jp*delta_u;

    if norm(pstar_Jp)<1e-6
        break;
    end
    %figure(count);
    [pstar_1_new,R1_new,g1_new] = findShape2(ustar(1:3*Ni),len1,Ni,R10,['Tube 1, iter ',num2str(count)]);
    [pstar_2_new,R2_new,g2_new] = findShape2(ustar(3*Ni+1:3*i*Ni),len1,Ni,R20,['Tube 2, iter ',num2str(count)]);
    %plotSection(g2_new,3e-3,[0 1 0]);
    %Compare updating with Jp and updating with findShape
    
    Rstar_2 = R2_new;
    Rstar_1 = R1_new;
    %pstar = [pstar_1_new; pstar_2_new];
    pstar = pstar_Jp+pstar;
    %disp(norm(pstar_Jp-pstar));



end

disp("Finished")
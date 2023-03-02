function res = boundary_condition(ini,tset)
    %ini is a initial guess which will be used for ode solving
    %tset is a constant tubeset type struct
    [u2out,a2out,~] = integrate_tube(ini,tset);
    
    u2z = u2out(3,:);
    a2 = a2out(:);
    res = zeros(2,1);
    
    alpha = tset.ctrl{2}(1) - tset.ctrl{1}(1);

    res(1) = (u2z(end) - 0);%namely uiz(L) = 0
    res(2) = (a2(1) - alpha); %namely a2(0) = theta2(0) - theta1(0)
    %res(3) = (u1z(end) - 0);
    %disp(res);
    %disp(res);
end
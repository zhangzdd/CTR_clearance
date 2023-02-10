function [pout,R_array,g_array]=findShape2(u,len,Ni,R0,label)
    %Find and plot shape of a Cosserat Model from stacked curvature values
    p = zeros(Ni,3);
    R_array = zeros(Ni,3,3);
    g_array = zeros(Ni,1,16);
    delta_s = len/Ni;
    R = R0;
    for m = 1:Ni     
        R_array(m,:,:) = R;
        p_dot = R*[0 0 1]';
        dR = R*hat(u(3*m-2:3*m))*delta_s;
        R = R+dR;
        %disp(R'*R);        
        if m == 1
            p(m,:) = 0;
        else
            p(m,:) = p(m-1,:)+p_dot'*delta_s;
        end
        g = [R p(m,:)';zeros(1,3) 1];
        g_array(m,:,:) = reshape(g,1,16);
    end

    plot3(p(:,1),p(:,2),p(:,3),"DisplayName",label);
    pout = [p(:,1)' p(:,2)' p(:,3)']';
%     disp(size(pout));
%     disp(size(p));
%     disp(reshape(pout,Ni,3)-p);
end
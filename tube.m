function t = tube(len,Ni,u,coef)
    t.len = len;% Length
    t.Ni = Ni;% Discretization value
    t.u = u;% Curvature function u(s)
    t.coef = coef;% Stiffness coefficient [kxy;poisson_ration]
    t.ds = t.len/t.Ni;% Length of a segment
    t.kxy = t.coef(1);
    t.kz = t.kxy/(t.coef(2)+1);
    t.K = diag([t.kxy, t.kxy, t.kz]);
end
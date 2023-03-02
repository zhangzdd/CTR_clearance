function coef = compute_coef_sum_V(c,Rstar_1,P2,Jp,S,pstar,Ni,circle_sample,lambda)
    coef = 0;
    for m = 1:circle_sample
        for n = 1:Ni
            k = (m-1)*Ni+n;
            
            Sj = squeeze(S(1,n,:,:));
            R1j = squeeze(Rstar_1(n,:,:));
            ci = c{m};
            P2j = squeeze(P2(n,:,:));

            coef = coef + lambda(k)*(ci'*R1j'*P2j*R1j*ci+2*ci'*R1j'*P2j*Sj*pstar+pstar'*Sj'*Sj*pstar);
        end
    end
end
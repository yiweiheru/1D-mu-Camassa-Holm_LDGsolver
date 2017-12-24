function [ Residue3 ] = Cresidue3(  x,Nelm,Ord,uh,rh,ph,Time )

elm_size=Ord+1;

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
un_der=zeros(elm_size,npt_quad);
uf=zeros(elm_size,2);

for k = 1 : npt_quad
    un(:,k)      =basis_1d(Ord,qpt(k));
    un_der(:,k)=basisDer_1d(Ord,qpt(k));
end

uf(:,1)=basis_1d(Ord,-1); %left boundary of elment
uf(:,2)=basis_1d(Ord,1);  %right boundary of element

Residue3=zeros(Nelm*elm_size,1);
Intergral=zeros(Nelm*elm_size,1);
Flux_L=zeros(Nelm*elm_size,1);
Flux_R=zeros(Nelm*elm_size,1);

%imply the periodic boundary condition to get the neighbourhood.
nei=zeros(Nelm,2);
for ne=1:Nelm
    if ne==1
        nei(ne,1)=Nelm;
        nei(ne,2)=2;
    elseif ne==Nelm
        nei(ne,1)=Nelm-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

%when implement the average of un, Jacobbi determinent is needed
mu_un=zeros(elm_size,1);
for ik=1:npt_quad
    mu_un(:,1)=mu_un(:,1)+qwt(ik)*un(:,ik);
end
mu_uh=0;
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for j=1:elm_size
        mu_uh=mu_uh+mu_un(j,1)*uh(ne,j)*Jaco;
    end
end

for ne=1:Nelm
    for ik=1:npt_quad
        u=0;p=0;r=0;
        for j=1:elm_size
            u=u+un(j,ik)*uh(ne,j);
            p=p+un(j,ik)*ph(ne,j);
            r=r+un(j,ik)*rh(ne,j);
        end
        f=2*mu_uh*u;
%         f=2*u;
        B=0.5*r^2;
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)+qwt(ik)*un_der(i,ik)*(-1*f+p-B);
        end
    end
 
    uLm=0;uLp=0;uRm=0;uRp=0;
    pLm=0;pRm=0;pLp=0;pRp=0;
    rLm=0;rRm=0;rLp=0;rRp=0;
    %compute flux
    for j=1:elm_size
        %f(uh)_hat take Lax_friedrichs flux
        uLm=uLm+uf(j,2)*uh(nei(ne,1),j);   %At (j-1/2), uh takes minus
        uLp=uLp+uf(j,1)*uh(ne,j);            %At (j-1/2), uh takes plus
        uRm=uRm+uf(j,2)*uh(ne,j);            %At (j+1/2), uh takes minus
        uRp=uRp+uf(j,1)*uh(nei(ne,2),j);   %At (j+1/2), uh takes plus
         
        %ph_hat takes minus
        pLm=pLm+uf(j,2)*ph(nei(ne,1),j); 
        pRm=pRm+uf(j,2)*ph(ne,j);
        pLp=pLp+uf(j,1)*ph(ne,j); 
        pRp=pRp+uf(j,1)*ph(nei(ne,2),j);
        
        %rh_hat takes minus
        rLm=rLm+uf(j,2)*rh(nei(ne,1),j);
        rRm=rRm+uf(j,2)*rh(ne,j);
        rLp=rLp+uf(j,1)*rh(ne,j);
        rRp=rRp+uf(j,1)*rh(nei(ne,2),j);
        
    end
        
    global orr
    orr = 1;%1:center; 2:LF
    if orr == 1
    fhat_L=0.5*(2*uLm+2*uLp)*mu_uh;
    fhat_R=0.5*(2*uRm+2*uRp)*mu_uh;
    elseif orr == 2
    fhat_L=0.5*(2*uLm+2*uLp-2*(uLp-uLm))*mu_uh;
    fhat_R=0.5*(2*uRm+2*uRp-2*(uRp-uRm))*mu_uh;  
    end
    
    phat_L=(pLm+pLp)/2;
    phat_R=(pRm+pRp)/2;
    
    Bhat_L=(0.5*rLm^2+0.5*rLp^2)/2;
    Bhat_R=(0.5*rRm^2+0.5*rRp^2)/2;

    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*(-1*fhat_L+phat_L-Bhat_L);
        Flux_R(num)=Flux_R(num)+uf(i,2)*(-1*fhat_R+phat_R-Bhat_R);
    end
end

Residue3=Residue3-(Intergral-Flux_R+Flux_L);

end



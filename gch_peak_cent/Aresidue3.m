function [ Residue3 ] = Aresidue3( x,Nelm,Ord,uh,rh,ph,Time )

elm_size=Ord+1;

npt_quad=Ord+3;
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
        B=0.5*r^2;
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)+qwt(ik)*un_der(i,ik)*(-1*f+p-B);
        end
    end
end


for ne=1:Nelm
%     if ne==1
%         uLp=1/26*(12*(0.5-Time)^2+23);
%         uLm=uLp;
%         uRp=uh(nei(ne,2),:)*uf(:,1);
%         uRm=uh(ne,:)*uf(:,2);
%         
%         rLm=12/13*(0.5-Time);
%         rRm=rh(ne,:)*uf(:,2);
%         
%         pLm=12/13*uLm+rLm*rLm;
%         pRm=ph(ne,:)*uf(:,2);
%     elseif ne==Nelm
%         uLp=uh(ne,:)*uf(:,1);
%         uLm=uh(nei(ne,1),:)*uf(:,2);
%         uRp=1/26*(12*(0.5-Time)^2+23);
%         uRm=uRp;
%         
%         rLm=rh(nei(ne,1),:)*uf(:,2);
%         rRm=12/13*(0.5-Time);
%         
%         pLm=ph(nei(ne,1),:)*uf(:,2);
%         pRm=12/13*uRm+rRm*rRm;
%     else
        uLp=uh(ne,:)*uf(:,1);                   
        uLm=uh(nei(ne,1),:)*uf(:,2);
        uRp=uh(nei(ne,2),:)*uf(:,1);          
        uRm=uh(ne,:)*uf(:,2);
        
        rLm=rh(nei(ne,1),:)*uf(:,2);
        rRm=rh(ne,:)*uf(:,2);
        
        pLm=ph(nei(ne,1),:)*uf(:,2);
        pRm=ph(ne,:)*uf(:,2);
%     end
    
    fhat_L=0.5*(2*(uLm+uLp)*mu_uh-2*(uLp-uLm)*abs(mu_uh));
    fhat_R=0.5*(2*(uRm+uRp)*mu_uh-2*(uRp-uRm)*abs(mu_uh));
    
    phat_L=pLm;
    phat_R=pRm;
    
    Bhat_L=0.5*rLm^2;
    Bhat_R=0.5*rRm^2;
    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*(-1*fhat_L+phat_L-Bhat_L);
        Flux_R(num)=Flux_R(num)+uf(i,2)*(-1*fhat_R+phat_R-Bhat_R);
    end
end

Residue3=Residue3-(Intergral-Flux_R+Flux_L);

end



function [ Residue2 ] = residue2( x,Nelm,Ord,uh,rh )

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

uf(:,1)=basis_1d(Ord,-1);
uf(:,2)=basis_1d(Ord,1);

% get the mean of uh
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


Residue2=zeros(Nelm*elm_size,1);
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

for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        r=0;
        for j=1:elm_size
            r=r+un(j,ik)*rh(ne,j);
        end
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)-qwt(ik)*un_der(i,ik)*r;
        end
    end
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Intergral(num)=Intergral(num)-Jaco*mu_un(i,1)*mu_uh;
    end
    
        
    rLm=0;rLp=0;rRm=0;rRp=0;
    for j=1:elm_size
        rLm=rLm+uf(j,2)*rh(nei(ne,1),j);
        rLp=rLp+uf(j,1)*rh(ne,j);
        rRm=rRm+uf(j,2)*rh(ne,j);
        rRp=rRp+uf(j,1)*rh(nei(ne,2),j);
    end
    
    rhat_L=0.5*(rLm+rLp);
    rhat_R=0.5*(rRm+rRp);
    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*rhat_L;
        Flux_R(num)=Flux_R(num)+uf(i,2)*rhat_R;
    end
end

Residue2=Residue2-(Intergral+Flux_R-Flux_L);

end



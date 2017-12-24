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
    for ik=1:npt_quad
        u=0;r=0;
        for j=1:elm_size
            u=u+un(j,ik)*uh(ne,j);
            r=r+un(j,ik)*rh(ne,j);
        end
        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num)+qwt(ik)*un_der(i,ik)*u*r;
        end
    end
     
    uLp=0;uRp=0;
    rLm=0;rLp=0;rRm=0;rRp=0;
    for j=1:elm_size
        uLp=uLp+uf(j,1)*uh(ne,j);                   %At (j-1/2), uh takes plus
        uRp=uRp+uf(j,1)*uh(nei(ne,2),j);    	    % At (j+1/2), uh takes plus
        
        rLm=rLm+uf(j,2)*rh(nei(ne,1),j);          %At (j-1/2), rh takes minus
        rLp=rLp+uf(j,1)*rh(ne,j);                     %At (j-1/2), rh takes plus
        rRm=rRm+uf(j,2)*rh(ne,j);                   %At (j+1/2), rh takes minus
        rRp=rRp+uf(j,1)*rh(nei(ne,2),j);            %At (j+1/2), rh takes plus
    end
    uhat_L=uLp;
    uhat_R=uRp;
    
    bhat_L=0.5*(rLm+rLp);
    bhat_R=0.5*(rRm+rRp);
    
    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*uhat_L*bhat_L;
        Flux_R(num)=Flux_R(num)+uf(i,2)*uhat_R*bhat_R;
    end
end

Residue2=Residue2-(Intergral-Flux_R+Flux_L);

end



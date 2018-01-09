function [ U_mod ] = Limiter( U,Ord,Nelm,x )
%LIMITER 此处显示有关此函数的摘要

elm_size=Ord+1;
uh=uhTransform(Nelm,elm_size,U);

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
for k = 1 : npt_quad
    un(:,k)      =basis_1d(Ord,qpt(k));
end

mu_un=zeros(elm_size,1);
for ik=1:npt_quad
    mu_un(:,1)=mu_un(:,1)+qwt(ik)*un(:,ik);
end
mu_uh=zeros(1,Nelm);
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for j=1:elm_size
        mu_uh(1,ne)=mu_uh(1,ne)+mu_un(j,1)*uh(ne,j)*Jaco;
    end
    mu_uh(1,ne)=mu_uh(1,ne)/(x(ne+1)-x(ne));
end


U_mod=U;
for ne=2:Nelm-1
    ur=uh(ne,elm_size)-mu_uh(1,ne);
    ul=mu_uh(1,ne)-uh(ne,1);
    
    ur_mod=minmod(ur,mu_uh(1,ne+1)-mu_uh(1,ne),mu_uh(1,ne)-mu_uh(1,ne-1));
    ul_mod=minmod(ul,mu_uh(1,ne+1)-mu_uh(1,ne),mu_uh(1,ne)-mu_uh(1,ne-1));
    
    U_mod(ne*elm_size)=mu_uh(1,ne)+ur_mod;
    U_mod((ne-1)*elm_size+1)=mu_uh(1,ne)-ul_mod;
end

end


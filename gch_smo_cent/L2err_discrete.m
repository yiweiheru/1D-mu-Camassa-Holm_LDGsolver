function [ l2error ] = L2err_discrete( Ustore,ir_num,ord_num,period )

NelmX=10*2^(ord_num-1);
dx=period/NelmX;
X=-period/2:dx:period/2;
num_pt=20;
xk=zeros(1,NelmX*(num_pt+1));

for ne=1:NelmX
    dx=X(ne+1)-X(ne);
    indent_bnd=dx/(num_pt+1)^2;
    temp=linspace(X(ne),X(ne+1),num_pt+1);
    xk((ne-1)*(num_pt+1)+2:ne*(num_pt+1)-1)=temp(2:num_pt);
    xk((ne-1)*(num_pt+1)+1)=temp(1)+indent_bnd;
    xk(ne*(num_pt+1))=temp(num_pt+1)-indent_bnd;
end

evalU=zeros(NelmX*(num_pt+1),ir_num,ord_num);

for Ord=1:ord_num
    for ir=1:ir_num
        Nelm=10*2^(ir-1);
        dx=period/Nelm;
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        sizeU=elm_size*Nelm;
        U=Ustore(1:sizeU,ir,Ord);
        uh=uhTransform(Nelm,Ord+1,U);
        for ne=1:NelmX
            for i=1:num_pt+1
                evalU((ne-1)*(num_pt+1)+i,ir,Ord)=...
                    evalue_uh( Ord,Nelm,x,uh,xk((ne-1)*(num_pt+1)+i) );
            end
        end
    end
end

l2error=zeros(ord_num,ir_num-1);
for Ord=1:ord_num
    for ir=1:ir_num-1
        dif=evalU(:,ir+1,Ord)-evalU(:,ir,Ord);
        val=0;
        NP=0;
        for n=1:NelmX*(num_pt+1)
                val=val+(dif(n))^2;
                NP=NP+1;
        end
        l2error(Ord,ir)=sqrt(val/NP);
    end
end

end


function [ l2error ] = LMerr_discrete( UStore,UexcStore,ir_num,ord_num,period )

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
evalUexc=zeros(NelmX*(num_pt+1),ir_num,ord_num);

for Ord=1:ord_num
    for ir=1:ir_num
        Nelm=10*2^(ir-1);
        dx=period/Nelm;
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        sizeU=elm_size*Nelm;
        U=UStore(1:sizeU,ir,Ord);
        Uexc=UexcStore(1:sizeU,ir,Ord);        
        uh=uhTransform(Nelm,Ord+1,U);
        uhexc=uhTransform(Nelm,Ord+1,Uexc);
        for ne=1:NelmX
            for i=1:num_pt+1
                evalU((ne-1)*(num_pt+1)+i,ir,Ord)=...
                    evalue_uh( Ord,Nelm,x,uh,xk((ne-1)*(num_pt+1)+i) );
                evalUexc((ne-1)*(num_pt+1)+i,ir,Ord)=...
                    evalue_uh( Ord,Nelm,x,uhexc,xk((ne-1)*(num_pt+1)+i) );                
            end
        end
    end
end

l2error=zeros(ord_num,ir_num);

for Ord=1:ord_num
    for ir=1:ir_num
        Du=evalU(:,ir,Ord)-evalUexc(:,ir,Ord);
        val=0;
        for n=1:NelmX*(num_pt+1)
            if (xk(n)<=-5 && xk(n)>=-25) || (xk(n)>=5 && xk(n)<=25)
                if abs(Du(n))>=val
                    val=abs(Du(n));
                end
            end
        end
        l2error(Ord,ir)=val;
    end
end

end

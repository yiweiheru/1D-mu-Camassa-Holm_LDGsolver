ord_num=2;
ir_num=6;
period=2*pi;

UStore=zeros((ord_num+1)*(10*2^(ir_num-1))+1,ir_num,ord_num);
UexcStore=zeros((ord_num+1)*(10*2^(ir_num-1))+1,ir_num,ord_num);
for Ord=1:ord_num
    for ir=1:ir_num
        Nelm=10*2^(ir-1);
        dx=period/Nelm;
        x=-period/2:dx:period/2;
        
        elm_size=Ord+1;
        P0=zeros(Nelm*elm_size,1);
        for ne=1:Nelm
            for i=1:elm_size
                xtemp=x(ne)+(x(ne+1)-x(ne))*(i-1)/(elm_size-1);
                num=(ne-1)*elm_size+i;
                P0(num,1)=2*sin(xtemp);
            end
        end
        
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        
        U=(massMat_inv*Amat)\P0;
        lengthU=elm_size*Nelm;
        UStore(1:lengthU,ir,Ord)=U;
        
        Uexc=P0/2;
        UexcStore(1:lengthU,ir,Ord)=Uexc;
    end
end
[ L2error ] = L2err_discrete( UStore,UexcStore,ir_num,ord_num,period );
for Ord=1:ord_num
    disp(ErrorOrder(L2error(Ord,:)));
end
% figure(1)
% plot_uh( U,Ord,Nelm,x,-period/2,period/2 )
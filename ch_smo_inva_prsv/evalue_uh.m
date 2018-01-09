function Val= evalue_uh( Ord,Nelm,x,uh,x0 )

for ne=1:Nelm
    if x0>=x(ne) && x0<x(ne+1)
        xi=(x0-x(ne))/(x(ne+1)-x(ne))*2+(-1);
        un=basis_1d(Ord,xi);
        Val=uh(ne,:)*un;
    elseif x0==x(Nelm+1)
        un=basis_1d(Ord,1);
        Val=uh(Nelm,:)*un;
    end
end

end
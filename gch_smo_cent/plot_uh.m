function plot_uh( U,Ord,Nelm,x ,XL,XR)

uh=uhTransform(Nelm,Ord+1,U);
xk=zeros(1,Nelm*(Ord+1));

for ne=1:Nelm
    dx=x(ne+1)-x(ne);
    indent_bnd=dx/(Ord+1)^2;
    temp=linspace(x(ne),x(ne+1),Ord+1);
    xk((ne-1)*(Ord+1)+2:ne*(Ord+1)-1)=temp(2:Ord);
    xk((ne-1)*(Ord+1)+1)=temp(1)+indent_bnd;
    xk(ne*(Ord+1))=temp(Ord+1)-indent_bnd;
end

it=0;
for np=1:Nelm*(Ord+1)
    if xk(np)>=XL && xk(np)<=XR
        it=it+1;
        x_to_plot(1,it)=xk(np);
        uh_to_plot(1,it)=evalue_uh(Ord,Nelm,x,uh,xk(np));
    end
end

if Ord==1 
plot(x_to_plot,uh_to_plot,'-r')
else
plot(x_to_plot,uh_to_plot,'-k')
end


end



function [ u ] = minmod( u1,u2,u3 )
%MINMOD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if sign(u1)==sign(u2) && sign(u2)==sign(u3)
    s=sign(u1);
    u=s*min([abs(u1),abs(u2),abs(u3)]);
else
    u=0;
end

end


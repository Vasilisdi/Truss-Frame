function [ D ] = Calc_det( A )
if size(A,1)~=size(A,2)
    message='Error:Matrix is not square';
    disp(message)
else 
    if max(size(A))==2
        D=A(1,1)*A(2,2)-A(1,2)*A(2,1);
    else
        for i=1:size(A,1)
            A_temp=A;
            A_temp(1,:)=[];
            A_temp(:,i)=[];
            if i==1
                D = (A(1,i)*(((-1)^(1+i))*Calc_det(A_temp)));
            else
                D=D+(A(1,i)*(((-1)^(1+i))*Calc_det(A_temp)));
            end
        end
    end
    
end
end
function B = Calc_inv(A)
% Dimensions of A
[rows,columns] = size(A);
% Keep on under square condition
if rows ~= columns
 disp('Only Square Matrices, please')
 B = [];
return
end
% Gauss-Jordan method
% Hence : [A|B_temp]-> [I|B_temp=inverse of A], determine B_temp=I =>
% =>inv(A)=B_temp
B_temp = eye(rows);
% Gauss-Jordan inversation
for j = 1 : rows
for i = j : rows
if A(i,j) ~= 0
for k = 1 : rows
        s = A(j,k); A(j,k) = A(i,k); A(i,k) = s;
        s = B_temp(j,k); B_temp(j,k) = B_temp(i,k); B_temp(i,k) = s;
end
        temp = 1/A(j,j);
for k = 1 : rows
        A(j,k) = temp * A(j,k);
        B_temp(j,k) = temp * B_temp(j,k);
end
for l = 1 : rows
if  l ~= j
        temp = -A(l,j);
for k = 1 : rows
        A(l,k) = A(l,k) + temp * A(j,k);
        B_temp(l,k) = B_temp(l,k) + temp * B_temp(j,k);
end
end
end
end
break
end
% Display warning if a row full of zeros is found
if A(i,j) == 0
disp('Warning: Singular Matrix')
B = 'error';
return
end
B=B_temp;
end
end
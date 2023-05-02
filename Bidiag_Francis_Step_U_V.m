function [U,T_Diagonal,V] = Bidiag_Francis_Step( U,B,V )
% T_spectrum will store the eigenvalues of T in its diagonal
T_Diagonal = B;
n = length(B);
% epsilon is the tolerance of convergence of the sub-diagonal element
eps=1e-18;
% in every iteration i we keep on applying francis steps until we get 0
% in the T_spectrum(i,i-1) entry, the i^th sub-diagonal element.
%For loop runs from n to 3 and not from n to 2 or n to 1.
for i=n:-1:3
% check for convergence of the i^th sub-diagonal entry
while T_Diagonal(i,i-1) > eps * sqrt(abs(T_Diagonal(i-1,i-1))+abs(T_Diagonal(i,i))) 
    %if not converged, apply francis step again.
    T_Diagonal(1:i,1:i) = Francis_Step(T_Diagonal(1:i,1:i)); 
    %modify this copy so that you also update UA and VA by applying the Givens  rotations, 
    % thus accumulating the U and V such that A = UΣVT
    U = Givens_rotation(U);
    V = Givens_rotation(V);
end
end
% In the base case we use the built-in eig function to compute the
% eigen values of the remaining 2x2 matrix in the UL corner.
T_Diagonal(1,2)=T_Diagonal(2,1); %This step is very important
eigs_of_2x2 = eig(T_Diagonal(1:2,1:2));
T_Diagonal(1,1) = eigs_of_2x2(1);
T_Diagonal(2,2) = eigs_of_2x2(2);
U = Givens_rotation(U);
V = Givens_rotation(V);
end
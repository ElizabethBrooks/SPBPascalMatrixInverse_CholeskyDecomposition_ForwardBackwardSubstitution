%{
  Name: Elizabeth Brooks
  File: SPDPascalMatrixInverse_CholeskyDecomposition_ForwardBackwardSubstitution
  Modified: 31 October 2017
%}
%Main function of script
function SPDPascalMatrixInverse_CholeskyDecomposition_ForwardBackwardSubstitution
%%Pascal matrix factorization
pn = 4 %For pn > 20, need to force symmetry using 1/2(B+B')
pi = 20.3 %Initialize scalar for Pascal matrix A(p,pn)
P = pascalMatrixA(pi,pn)
C = choleskyA(P,pn)
%C*C' %Use to verify results
det = determinantA(C,pn)
%%Pascal matrix Inverse
in = pn
N = inv(A)
N = computeInverseA(C,in)
CN = choleskyA(N,in) %Check if SPD using cholesky
detN = determinantA(CN,pn)
%
%Compute pascal matrix for integer n and scalar p
% at O(n^2) efficiency
function A = pascalMatrixA(p,n)
  A = zeros(n); %Initialize matrix A with zeros
  for j = 1:n %Columnwise on A
    for i = 1:n %Rowwise on A
      if((i == 1) || (j == 1))
        A(i,j) = p; %Set the first column or row to p
      else
        A(i,j) = A(i,j-1)+A(i-1,j); %Set remaining values
      end
    end
  end
end %End function pascalMatrixA
%
%Cholesky factorization A = LL' at O((n^3)/3) efficiency
function A = choleskyA(A,n)
for j = 1:n %Columnwise on A
  for i = 1:n %Rowwise on A
    if(i < j) %Above diagonal
      A(i,j) = 0;
    elseif(i == j) %Diagonal
      if(A(i,j) <= 0) %Check if SPD
        Err = 'Error: not SPD matrix.' %Display error message
        return; %Exit matrix computations
      end
      A(i,j) = sqrt(A(i,j)); %Corner
    else %Below diagonal
      A(i,j) = A(i,j)/A(j,j); %l
    end
  end
  for z = j+1:n %Columnwise on submatrix of A
    for i = z:n %Rowwise on submatrix of A
      A(i,z) = A(i,z) - (A(i,j)*A(z,j)); %a = ll'
    end
  end
end
end %End function choleskyA
%
%Compute determinant of square SPD matrix
% as the square of the products of the diagonal elements of L
function det = determinantA(A,n)
det = 1; %Initialize determinant
for j = 1:n %Index on L
  det = det*A(j,j); %The products of the diagonal elements of L
end
det = det*det; %Square the products of the diagonal elements of L
end %End function determinantA
%
%Solve systems of linear equations with forward and backward substitution
%Forward substitution
function y = forwardSubstitutionA(A,b,n)
y = zeros(n,1); %Matrix to be solved, y
for j = 1:n %Columnwise on L
  y(j) = b(j)/A(j,j); %Matrix is SPD, no need to check A(j,j) > 0
  for i = j+1:n %Rowwise on L
    if(A(i,j) ~= 0) %Skip computations with zero multiplication
      b(i) = b(i) - (A(i,j)*y(j));
    end
  end
end
end %End function forwardSubstitutionA
%
%Backward substitution
function x = backwardSubstitutionA(A,y,n)
x = zeros(n,1); %Matrix to be solved, x
for j = n:-1:1 %Decrement columnwise index on L'
  x(j) = y(j)/A(j,j); %Matrix is SPD, no need to check A(j,j) > 0
  for i = 1:j-1 %Rowwise on L'
    if(A(j,i) ~= 0) %Skip computations with zero multiplication
      y(i) = y(i) - (A(j,i)*x(j));
    end
  end
end
end %End function backwardSubstitutionA
%
%Compute the inverse of input matrix A of input size n
function N = computeInverseA(A,n)
  b = zeros(n,n); %Input matrix, b
  y = zeros(n,1); %Initialize vector y
  for j =1:n
    for i = 1:n
      if(i == j)
        b(i,j) = 1;
      end
    end
  end
  for j = 1:n
    for i = 1:n
      y = forwardSubstitutionA(A,b(i,:),n);
      N(i,:) = backwardSubstitutionA(A,y,n);
    end
  end
end %End function computeInverseA
%
end %End main function of script
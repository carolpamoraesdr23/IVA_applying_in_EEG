function [h,invQ,R]=decouple_trick(W,n,invQ,R)
% h=decouple_trick(W,n)
% h=decouple_trick(W,n,1)
% [h,invQ]=decouple_trick(W,n,invQ)
% [h,Q,R]=decouple_trick(W,n,Q,R)
%
% Computes the h vector for the decoupling trick [1] of the nth row of W. W
% can be K 'stacked' square matrices, i.e., W has dimensions N x N x K.
% The output vector h will be formatted as an N x K matrix.  There are many
% possible methods for computing h.  This routine provides four different
% (but of course related) methods depending on the arguments used.
%
% Method 1:
% h=decouple_trick(W,n)
% h=decouple_trick(W,n,0)
% -Both calls above will result in the same algorithm, namely the QR
% algorithm is used to compute h.
%
% Method 2:
% h=decouple_trick(W,n,~), where ~ is anything
% -Calls the projection method.
%
% Method 3:
% [h,invQ]=decouple_trick(W,n,invQ)
% -If two output arguments are specified then the recursive algorithm
% described in [2].  It is assumed that the decoupling will be performed in
% sequence see the demo subfunction for details.
% An example call sequence:
%  [h1,invQ]=decouple_trick(W,1);
%  [h2,invQ]=decouple_trick(W,2,invQ);
%
% Method 4:
% [h,Q,R]=decouple_trick(W,n,Q,R)
% -If three output arguments are specified then a recursive QR algorithm is
% used to compute h.
% An example call sequence:
%  [h1,Q,R]=decouple_trick(W,1);
%  [h2,Q,R]=decouple_trick(W,2,Q,R);
%
% See the subfunction demo_decoupling_trick for more examples.  The demo
% can be executed by calling decouple_trick with no arguments, provides a
% way to compare the speed and determine the accuracy of all four
% approaches.
%
% Note that methods 2 & 3 do not normalize h to be a unit vector.  For
% optimization this is usually not of interest.  If it is then set the
% variable boolNormalize to true.
%
% Main References:
% [1] X.-L. Li & X.-D. Zhang, "Nonorthogonal Joint Diagonalization Free of Degenerate Solution," IEEE Trans. Signal Process., 2007, 55, 1803-1814
% [2] X.-L. Li & T. Adali, "Independent component analysis by entropy bound minimization," IEEE Trans. Signal Process., 2010, 58, 5151-5164
%
% Coded by Matthew Anderson (matt dot anderson at umbc dot edu)

% Version 01 - 20120919 - Initial publication


if nargin==0
   help decouple_trick
   demo_decouple_trick
   return
end
if nargin==1
   help decouple_trick
   error('Not enough inputs -- see displayed help above.')
end
[M,N,K]=size(W);
if M~=N
   error('Assuming W is square matrix.')
end
h=zeros(N,K);

% enables an additional computation that is usually not necessary if the
% derivative is of  interest, it is only necessary so that sqrt(det(W*W'))
% = sqrt(det(Wtilde*Wtilde'))*abs(w'*h) holds.  Furthermore, it is only
% necessary when using the recursive or projection methods.
%
% a user might wish to enable the calculation by setting the quantity below
% to true
boolNormalize=false;

if nargout==3
   % use QR recursive method
   % [h,Qnew,Rnew]=decouple_trick(W,n,Qold,Rold)
   if n==1
      invQ=zeros(N,N,K);
      R=zeros(N,N-1,K);
   end
   for k=1:K
      if n==1
         Wtilde=W(2:N,:,k);
         [invQ(:,:,k),R(:,:,k)]=qr(Wtilde');
      else
         n_last=n-1;
         e_last = zeros(N-1,1);
         e_last(n_last) = 1;
         [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),-W(n,:,k)',e_last);
         [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),W(n_last,:,k)',e_last);
      end
      h(:,k)=invQ(:,end,k); % h should be orthogonal to W(nout,:,k)'
   end
elseif nargout==2
   % use recursive method
   % [h,invQ]=decouple_trick(W,n,invQ), for any value of n=1, ..., N
   % [h,invQ]=decouple_trick(W,1), when n=1
   
   if n==1
      invQ=zeros(N-1,N-1,K);
   end
   % Implement a faster approach to calculating h.
   for k=1:K
      if n==1
         Wtilde=W(2:N,:,k);
         invQ(:,:,k)=inv(Wtilde*Wtilde');
      else
         if nargin<3
            help decouple_trick
            error('Need to supply invQ for recursive approach.')
         end
         [Mq,Nq,Kq]=size(invQ);
         if Mq~=(N-1) || Nq~=(N-1) || Kq~=K
            help decouple_trick
            error('Input invQ does not have the expected dimensions.')
         end
         n_last=n-1;
         Wtilde_last=W([(1:n_last-1) (n_last+1:N)],:,k);
         w_last=W(n_last,:,k)';
         w_current=W(n,:,k)';
         c = Wtilde_last*(w_last - w_current);
         c(n_last) = 0.5*( w_last'*w_last - w_current'*w_current );
         %e_last = zeros(N-1,1);
         %e_last(n_last) = 1;
         temp1 = invQ(:,:,k)*c;
         temp2 = invQ(:,n_last,k);
         inv_Q_plus = invQ(:,:,k) - temp1*temp2'/(1+temp1(n_last));
         
         temp1 = inv_Q_plus'*c;
         temp2 = inv_Q_plus(:,n_last);
         invQ(:,:,k) = inv_Q_plus - temp2*temp1'/(1+c'*temp2);
         % inv_Q is Hermitian
         invQ(:,:,k) = (invQ(:,:,k)+invQ(:,:,k)')/2;
      end
      
      temp1 = randn(N, 1);
      Wtilde = W([(1:n-1) (n+1:N)],:,k);
      h(:,k) = temp1 - Wtilde'*invQ(:,:,k)*Wtilde*temp1;
   end
   if boolNormalize
      h=vecnorm(h);
   end
elseif nargin==2 || invQ==0
   % use (default) QR approach
   % h=decouple_trick(W,n)
   % h=decouple_trick(W,n,0)
   for k=1:K
      [Q,zois]=qr(W([(1:n-1) (n+1:N)],:,k)');
      h(:,k)=Q(:,end); % h should be orthogonal to W(nout,:,k)'
   end
else % use projection method
   % h=decouple_trick(W,n,~), ~ is anything
   for k=1:K
      temp1 = randn(N, 1);
      Wtilde = W([(1:n-1) (n+1:N)],:,k);
      h(:,k) = temp1 - Wtilde'*((Wtilde*Wtilde')\Wtilde)*temp1;
   end
   if boolNormalize
      h=vecnorm(h);
   end
end

% function hessian_mat = hessian_sparse(func,x,varargin)
% Computes second order partial derivatives
% used the 7 points formula
%
% INPUTS
%    func:           name of the function
%    x:              vector of variables around which the Hessian is calculated
%    varargin:       list of arguments following x
%
% OUTPUTS
%    hessian_matrix: Hessian matrix
%
% ALGORITHM
% Uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
%
% SPECIAL REQUIREMENTS
%    none
%  
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.

%% computes second order partial derivatives
% uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
function hessian_mat = hessian_sparse(func,x,varargin)
global options_func = str2func(func);n=size(x,1);
%h1=max(abs(x),options_.gstep*ones(n,1))*eps^(1/3);h1=max(abs(x),sqrt(options_.gstep)*ones(n,1))*eps^(1/6);h_1=h1;xh1=x+h1;h1=xh1-x;xh1=x-h_1;h_1=x-xh1;xh1=x;f0=feval(func,x,varargin{:});nf = size(f0,1);f1=zeros(nf,n);f_1=f1;for i=1:n    xh1(i)=x(i)+h1(i);    f1(:,i)=feval(func,xh1,varargin{:});    xh1(i)=x(i)-h_1(i);    f_1(:,i)=feval(func,xh1,varargin{:});    xh1(i)=x(i);    i=i+1;endxh_1=xh1;hessian_mat = spalloc(nf,n*n,3*nf*n);for i=1:n    if i > 1        k=[i:n:n*(i-1)];        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);    end     hessian_mat(:,(i-1)*n+i)=sparse((f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i)));    temp=f1+f_1-f0*ones(1,n);        for j=i+1:n        xh1(i)=x(i)+h1(i);        xh1(j)=x(j)+h_1(j);        xh_1(i)=x(i)-h1(i);        xh_1(j)=x(j)-h_1(j);        hessian_mat(:,(i-1)*n+j)=sparse(-(-feval(func,xh1,varargin{:})-feval(func,xh_1,varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j)));        xh1(i)=x(i);        xh1(j)=x(j);        xh_1(i)=x(i);        xh_1(j)=x(j);        j=j+1;    end    i=i+1;end% 10/03/02 MJ used the 7 points formula
function [X,E,iter] = lrsr(A,B,lambda1,lambda2,max_iter, st_mat)

% Solve the Low-Rank and Sparse Representation (LRSR) minimization problem by M-ADMM
%
% min_{X,E} ||X||_*+lambda1*||X||_1+lambda2*loss(E), s.t. A=BX+E
% loss(E) = ||E||_1 or 0.5*||E||_F^2 or ||E||_{2,1}
% ---------------------------------------------
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

tol = 1e-8; 
rho = 1.2;
mu = 1e-4;
max_mu = 100;
DEBUG = 0;
loss = 'l21';

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'loss');        loss = opts.loss;            end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


[d,na] = size(A);
[~,nb] = size(B);

X = zeros(nb,na);
E = zeros(d,na);
Z = X;
J = X;

Y1 = E;
Y2 = X;
Y3 = X;
BtB = B'*B;
BtA = B'*A;
I = eye(nb);
invBtBI = (BtB+2*I)\I;
iter = 0;
for iter = 1 : max_iter
    Xk = X;
    Zk = Z;
    Ek = E;
    Jk = J;
    %%
    % first super block {Z,J,E}
    Z =  (mu*X+Y2) ./ (2*(st_mat .* st_mat) + mu); 
    %[Z,nuclearnormZ] = prox_nuclear(X+Y2/mu,1/mu);
    %% 
    %J = prox_l1(X+Y3/mu,lambda1/mu);
    J = prox_l1(X+Y3/mu,(lambda1/mu).*st_mat);
%     sum(sum(J))
%     J(J<0)=0;
    
    if strcmp(loss,'l1')
        E = prox_l1(A-B*X+Y1/mu,lambda2/mu);
    elseif strcmp(loss,'l21')
        E = prox_l21(A-B*X+Y1/mu,lambda2/mu);
    elseif strcmp(loss,'l2')
        E = mu*(A-B*X+Y1/mu)/(lambda2+mu);
    else
        error('not supported loss function');
    end
    % second  super block {X}
    X = invBtBI*(B'*(Y1/mu-E)+BtA-(Y2+Y3)/mu+Z+J);
    
  
    dY1 = A-B*X-E;
    dY2 = X-Z;
    dY3 = X-J;
    chgX = max(max(abs(Xk-X)));
    chgE = max(max(abs(Ek-E)));
    chgZ = max(max(abs(Zk-Z)));
    chgJ = max(max(abs(Jk-J)));
    chg = max([chgX chgE chgZ chgJ max(abs(dY1(:))) max(abs(dY2(:))) max(abs(dY3(:)))]);
    
    if chg < tol
        break;
    end 
    Y1 = Y1 + mu*dY1;
    Y2 = Y2 + mu*dY2;
    mu = min(rho*mu,max_mu);   
    
    
    if rem(iter,5) == 0
        disp(['Completed the ',num2str(iter), '/', num2str(max_iter), ' iteration of solution...'])
    end
end

%% 

function out = comp_loss(E,normtype)

switch normtype
    case 'l1'
        out = norm(E(:),1);
    case 'l21'
        out = 0;
        for i = 1 : size(E,2)
            out = out + norm(E(:,i));
        end
    case 'l2'
        out = 0.5*norm(E,'fro')^2;
end


function x = prox_l1(b,lambda)

% The proximal operator of the l1 norm
% 
% min_x lambda*||x||_1+0.5*||x-b||_2^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
x = max(0,b-lambda)+min(0,b+lambda);


function X = prox_l21(B,lambda)

% The proximal operator of the l21 norm of a matrix
% l21 norm is the sum of the l2 norm of all columns of a matrix 
%
% min_X lambda*||X||_{2,1}+0.5*||X-B||_2^2
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
%

X = zeros(size(B));
for i = 1 : size(X,2)
    nxi = norm(B(:,i));
    if nxi > lambda  
        X(:,i) = (1-lambda/nxi)*B(:,i);
    end
end

function [U,S,V] = svdecon(X)
% Input:
% X : m x n matrix
%
% Output:
% X = U*S*V'
%
% Description:
% Does equivalent to svd(X,'econ') but faster
%
[m,n] = size(X);
if  m <= n
    C = X*X';
    [U,D] = eig(C);
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    U = U(:,ix);    
    
    if nargout > 2
        V = X'*U;
        s = sqrt(d);
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    C = X'*X; 
    [V,D] = eig(C);
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    V = V(:,ix);    
    
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(d);
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end

clc;
clear;
close all;

%% Receive sequence
% homework
% m = 6;
% t = 4;
% r_coeff = [36 4 33 21 56 52 47 13 39 0 5 11 37];
% r_power = [0 1 2 3 4 5 6 7 8 9 10 23 60];

% textbook example 
m = 4;
t = 3;
r_coeff = [1 5 11];
r_power = [7 5 2];

%% Syndrome
S_terms = -Inf(2*t+1,length(r_power));
S = -Inf(2*t+1,1);
for i = 1:2*t
    S_terms(i,:) = i*r_power+r_coeff;
    S(i) = GF_add(S_terms(i,:),m);
end

%% Set up
s_r = -Inf(2*t,2*t+1);
t_r = -Inf(2*t,2*t+1);
Q = -Inf(2*t,2*t+1);
A = -Inf(4,2*t+1,2*t);
%% Initialization
s_r(1,2*t+1) = 0;
t_r(1,:) = S;
A([1,4],1,1) = 0;

%% Iterative procedure
for r = 1:2*t
    % Check degree
    deg_s = find(isfinite(s_r(r,:)),1,'last')-1;
    deg_t = find(isfinite(t_r(r,:)),1,'last')-1;
    deg_A = find(any(isfinite(A(:,:,r))),1,'last')-1;
    
    % Check recursive degree
    if deg_t < t
        break
    end
    
    % Recursive equation
    deg_Q = deg_s-deg_t;
    re = s_r(r,:);
    deg_re = deg_s;
    for i = deg_Q+1:-1:1
        Q(r,i) = mod((-1)*t_r(r,deg_t+1)+re(deg_re+1),2^m-1);
        for j = deg_t:-1:0
            % re = re-Q(i)*t
            re(i+j) = GF_add([re(i+j),Q(r,i)+t_r(r,j+1)],m);
            if re(deg_re+1) == -Inf
                deg_re = deg_re-1;
            end
        end     
    end
    
    % Calculate A^(r+1)
    A([1 2],:,r+1) = A([3 4],:,r);
    for i = 0:deg_A
        for j = 0:deg_Q
            A(3,i+j+1,r+1) = GF_add([A(3,i+j+1,r+1),Q(r,j+1)+A(3,i+1,r)],m);
            A(4,i+j+1,r+1) = GF_add([A(4,i+j+1,r+1),Q(r,j+1)+A(4,i+1,r)],m);
        end
    end
    for i = 1:2*t+1
        A(3,i,r+1) = GF_add([A(3,i,r+1),A(1,i,r)],m);
        A(4,i,r+1) = GF_add([A(4,i,r+1),A(2,i,r)],m);
    end
    
    % Calcualte t^(r+1)
    for i = 0:deg_t
        for j = 0:deg_Q
            t_r(r+1,i+j+1) = GF_add([t_r(r+1,i+j+1),Q(r,j+1)+t_r(r,i+1)],m);
        end
    end
    for i = 1:deg_s+1
        t_r(r+1,i) = GF_add([t_r(r+1,i),s_r(r,i)],m);
    end
    
    % Calculate s^(r+1)
    s_r(r+1,:) = t_r(r,:);
end

%% Calculate dr,Lambda(X) and Omega(X)
Lambda = -Inf(2*t,2*t+1);
for i = 2:r
    dr = A(4,1,i);
    dr_r = mod(-dr,2^m-1);
    Lambda(i,:) = GF_mul(dr_r+A(4,:,i),m);
    Omega = GF_mul(dr_r+t_r(i,:),m);
    deg_L = find(~isinf(Lambda(i,:)),1,'last')-1;

end

%% Error location
error_count = 0;
for i = 0:2^m-2
    Lambda_sum = -Inf;
    for j = 0:deg_L
        Lambda_sum = GF_add([Lambda_sum,Lambda(r,j+1)+i*j],m);
    end
    if Lambda_sum == -Inf
        error_count = error_count+1;
        Xi_r(error_count) = i;
    end
end
Xi_r_shift = circshift(Xi_r,[1,-1]);
Xi = mod(-Xi_r_shift,2^m-1);

%% Error evaluator <Omega>
Omega = -Inf(1,2*t);
for i = 1:2*t
    sum = Lambda(r,i)+S(1:2*t)';
    Omega = Inf_sum(Omega,[-Inf(1,i-1),sum(1:2*t+1-i)],m);
end

%% Error magnitude <Yi>
Omega_r = -Inf(1,deg_L);
A_mul = zeros(deg_L);
for i = 1:deg_L
    D = 0;
    Omega_r(i) = GF_add(Omega(1:deg_L)+Xi_r_shift(i)*(0:deg_L-1),m);
    
    for j = 1:deg_L
        if i ~= j
            A_mul(i,j) = GF_add([0,Xi_r_shift(i)+Xi(j)],m);
        else
            A_mul(i,j) = 0;
        end
        D = D+A_mul(i,j);
        DD = Xi(i)+D;
        Yi(i) = mod(Omega_r(i)+(-1)*DD,2^m-1);
    end       
end

%% Decode codeword
r_codeword = -Inf(1,2^m-1);
r_codeword(r_power+1) = r_coeff;
e_codeword = -Inf(1,2^m-1);
e_codeword(Xi+1) = Yi;
d_codeword = Inf_sum(r_codeword,e_codeword,m);
d_codeword(~isfinite(d_codeword)) = 0;

%% Print to .txt
fileID = fopen('q2.txt','w');
fprintf(fileID,'Px[a0,a1,a2,...ai,..] denotes the polynomial martix, where i is the power of X and ai is the power of the primitive element.\n');
fprintf(fileID,'If ai= -inf,then (alpha)^(ai)= 0 where alpha is primitive element in GF(2^m). And if ai= 0,then (alpha)^(ai)= 1.\n');
fprintf(fileID,'\n');
for i = 1:deg_L
    fprintf(fileID,'-----------Euclidean_algorithm Iteration %d-----------\n',i);
    % r
    fprintf(fileID,'r = %d\n',i);
    % s(r) (x)
    format = ['s(r) (x) = Px[' repmat(' %1.0f ',1,length(s_r(i+1,~isinf(s_r(i+1,:))))) ']\n'];
    fprintf(fileID,format,s_r(i+1,~isinf(s_r(i+1,:))));
    % t(r) (x)
    format = ['t(r) (x) = Px[' repmat(' %1.0f ',1,length(t_r(i+1,~isinf(t_r(i+1,:))))) ']\n'];
    fprintf(fileID,format,t_r(i+1,~isinf(t_r(i+1,:))));
    % Q(r) (x)
    format = ['Q(r) (x) = Px[' repmat(' %1.0f ',1,length(Q(i,~isinf(Q(i,:))))) ']\n'];
    fprintf(fileID,format,Q(i,~isinf(Q(i,:))));
    % Lambda
    format = ['Lambda (x) = Px[' repmat(' %1.0f ',1,length(Lambda(i+1,~isinf(Lambda(i+1,:))))) ']\n'];
    fprintf(fileID,format,Lambda(i+1,~isinf(Lambda(i+1,:))));
    % A(r) (x)
    fprintf(fileID,'A(r)(x) = [A11(x) A12(x) ; A21(x) A22(x)]\n');
    format = ['A11(r) (x) = Px[' repmat(' -Inf ',1,length(A(2,~isinf(A(2,:,i+1)),i+1))) ']\n'];
    fprintf(fileID,format);
    format = ['A12(r) (x) = Px[' repmat(' %1.0f ',1,length(A(2,~isinf(A(2,:,i+1)),i+1))) ']\n'];
    fprintf(fileID,format,A(2,~isinf(A(2,:,i+1)),i+1));
    format = ['A21(r) (x) = Px[' repmat(' -Inf ',1,length(A(2,~isinf(A(4,:,i+1)),i+1))) ']\n'];
    fprintf(fileID,format);
    format = ['A22(r) (x) = Px[' repmat(' %1.0f ',1,length(A(2,~isinf(A(4,:,i+1)),i+1))) ']\n'];
    fprintf(fileID,format,A(4,~isinf(A(4,:,i+1)),i+1));
    fprintf(fileID,'\n');
    
end
X_format = ['Error location Xi = ' repmat(' X^%1.0f ',1,length(Xi)) '\n'];
fprintf(fileID,X_format,Xi);
Y_format = ['Error magnitude Yi = ' repmat(' alpha^%1.0f ',1,length(Xi)) '\n'];
fprintf(fileID,Y_format,Yi);
fprintf(fileID,'\n');
d_format = ['Decoded codeword = (' repmat('%1.0f ',1,length(d_codeword)) ')\n'];
fprintf(fileID,d_format,d_codeword);
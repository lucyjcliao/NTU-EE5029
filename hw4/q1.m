clc;
clear;
close all;

%% Receive sequence
% homeword
m = 6;
t = 4;
r_coeff = [36 4 33 21 56 52 47 13 39 0 5 11 37];
r_power = [0 1 2 3 4 5 6 7 8 9 10 23 60];

% textbook example
% m = 4;
% t = 3;
% r_coeff = [1 5 11];
% r_power = [7 5 2];

%% Syndrome
S_terms = -Inf(2*t,length(r_power));
S = -Inf(2*t,1);
for i = 1:2*t
    S_terms(i,:) = i*r_power+r_coeff;
    S(i) = GF_add(S_terms(i,:),m);
end

%% Set up
dr = zeros(2*t+1,1);
T = -Inf(2*t+1,2*t);
B = -Inf(2*t+1,2*t);
Lambda = -Inf(2*t+1,2*t);
L = zeros(2*t+1,1);
%% Initialization
dr(1) = 0;
B(1,1) = 0;
Lambda(1,1) = 0;

%% Iterative procedure
for r = 1:2*t
    % Calculate dr
    dr_terms = GF_mul([Lambda(r,:);circshift(flip(S),r)'],m);
    dr(r+1) = GF_add(dr_terms(~isinf(dr_terms)),m);
    
    % Check <dr==0>
    if isinf(dr(r+1))
        T(r+1,:) = Lambda(r,:);         
    else
        B_shift = circshift(B(r,:),[1,1]);
        for i = 1:2*t
            T(r+1,i) = GF_add([Lambda(r,i);dr(r+1)+B_shift(i)],m);
        end
    end
    
    % Check <2*L <= r-1>
    if 2*L(r) <= r-1
        B(r+1,:) = mod(-dr(r+1),2^m-1)+Lambda(r,:);
        for j = 1:2*t
            B(r+1,j) = GF_add(B(r+1,j),m);
        end
        Lambda(r+1,:) = T(r+1,:);
        L(r+1) = r-L(r);
    else
        B(r+1,:) = circshift(B(r,:),[1,1]);
        Lambda(r+1,:) = T(r+1,:);
        L(r+1) = L(r);
    end
end

%% Check <deg(Lambda(X)) = L>
deg_Lambda = sum(~isinf(Lambda(r+1,:)))-1;
if deg_Lambda > L(r+1)
    error('More than t errors');
end

%% Error locations <Xi>
error_count = 0;
for i = 0:2^m-2
    Lambda_sum = -Inf;
    for j = 0:L(r+1)
        Lambda_sum = GF_add([Lambda_sum,Lambda(r+1,j+1)+i*j],m);
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
    sum = Lambda(r+1,i)+S(1:2*t)';
    Omega = Inf_sum(Omega,[-Inf(1,i-1),sum(1:2*t+1-i)],m);
end


%% Error magnitude <Yi>
Omega_r = -Inf(1,L(r+1));
D = zeros(1,L(r+1));
for i = 1:L(r+1)
    Omega_r(i) = GF_add(Omega(1:L(r+1))+Xi_r_shift(i)*(0:L(r+1)-1),m);
    
    for j = 1:L(r+1)
        if i ~= j
            Xi_r_shift(i);
            Xi(j);
            A(i,j) = GF_add([0,Xi_r_shift(i)+Xi(j)],m);
        else
            A(i,j) = 0;
        end
        D(i) = D(i)+A(i,j);
        DD(i) = Xi(i)+D(i);
        Yi(i) = mod(Omega_r(i)+(-1)*DD(i),2^m-1);
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
fileID = fopen('q1.txt','w');
fprintf(fileID,'Px[a0,a1,a2,...ai,..] denotes the polynomial martix, where i is the power of X and ai is the power of the primitive element.\n');
fprintf(fileID,'If ai= -inf,then (alpha)^(ai)= 0 where alpha is primitive element in GF(2^m). And if ai= 0,then (alpha)^(ai)= 1.\n');
fprintf(fileID,'\n');
for i = 1:2*t
    fprintf(fileID,'-----------BM_algorithm Iteration %d-----------\n',i);
    % r
    fprintf(fileID,'r = %d\t',i);
    % dr
    fprintf(fileID,'dr = %d\n',dr(i+1));

    format = ['Px[' repmat(' %1.0f ',1,2*t) ']\n'];
    % Lambda(X)
    fprintf(fileID,'Lambda(X) = ');
    fprintf(fileID,format,Lambda(i+1,:));
    % B(X)
    fprintf(fileID,'B(X) = ');
    fprintf(fileID,format,B(i+1,:));
    % T(X)
    fprintf(fileID,'T(X) = ');
    fprintf(fileID,format,T(i+1,:));
    % L
    fprintf(fileID,'L = %d\n',L(i+1));
    fprintf(fileID,'\n');
end
X_format = ['Error location Xi = ' repmat(' X^%1.0f ',1,length(Xi)) '\n'];
fprintf(fileID,X_format,Xi);
Y_format = ['Error magnitude Yi = ' repmat(' alpha^%1.0f ',1,length(Xi)) '\n'];
fprintf(fileID,Y_format,Yi);
fprintf(fileID,'\n');
d_format = ['Decoded codeword = (' repmat('%1.0f ',1,length(d_codeword)) ')\n'];
fprintf(fileID,d_format,d_codeword);


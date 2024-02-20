%% experimental setup
% Noiseless Case
%% N = 20
N = 20;
M = N;
S_max = N;
success_prob_20 = zeros(S_max,M);
avg_error_20 = zeros(S_max,M);
for m = 1:N
    for s = 1:S_max
        success = 0;
        err = 0;
        for i = 1:2000
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
    
            x = zeros(1,N); % generate sparse vector 
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';
            y = A*x';
            x_approx = OMP(y,A,N);
        
             if(vecnorm(x_approx-x)/vecnorm(x) <= 10^-5)
                 success = success+1;
             else
                 err = err+1;
             end
        end
        success_prob_20(s,m) = success;
        avg_error_20(s,m) = err;
    end
end
%% generate plots
figure(1)
imagesc(success_prob_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 20']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(2)
imagesc(avg_error_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 20']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar
%% N = 50
N = 50;
M = N;
S_max = N;
success_prob_50 = zeros(S_max,M);
avg_error_50 = zeros(S_max,M);
for m = 1:N
    for s = 1:S_max
        success = 0;
        err = 0;
        for i = 1:2000
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
    
            x = zeros(1,N); % generate sparse vector 
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';
            y = A*x';
            x_approx = OMP(y,A,N);
        
             if(vecnorm(x_approx-x)/vecnorm(x) <= 10^-5)
                 success = success+1;
             else
                 err = err+1;
             end
        end
        success_prob(s,m) = success;
        avg_error(s,m) = err;
    end
end
%%
figure(3)
imagesc(success_prob_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 50']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(4)
imagesc(avg_error_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 50']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

%% N = 100
N = 100;
M = N;
S_max = N;
success_prob_100 = zeros(S_max,M);
avg_error_100 = zeros(S_max,M);
for m = 1:N
    for s = 1:S_max
        success = 0;
        err = 0;
        for i = 1:2000
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
    
            x = zeros(1,N); % generate sparse vector 
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';
            y = A*x';
            x_approx = OMP(y,A,N);
        
             if(vecnorm(x_approx-x)/vecnorm(x) <= 10^-5)
                 success = success+1;
             else
                 err = err+1;
             end
        end
        success_prob_100(s,m) = success;
        avg_error_100(s,m) = err;
    end
end
%%
figure(5)
imagesc(success_prob_100/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 100']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(6)
imagesc(avg_error_100/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 100']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar
%%
% OMP
% input signal - y; residual - r; approximation - a_x; support vector - Lambda
% initial setup: r_0 = s, a_0 = 0,

% The OMP function takes input signal y, measurement matrix A, and the
% stopping threshold thrs
function a_x = OMP(y, A, thrs)
    % initial setup
    [~,n] = size(A);
    %Lambda = zeros(thrs,thrs);
    r = y;
    A_new = [];
    a_x = zeros(1,n);
    i = 1;
    while (i <= thrs)
        % finding index of column of the matrix A that has the maximum
        % projection from the residual 
        [~,colm] = max(abs(A'*r));
      
        % updating the list of atoms 
        A_new = [A_new A(:,colm)];
        
        % asign 0 to the used column vector
        A(:,colm) = 0;

        % obtain lambda using least square
        lambda = pinv(A_new)*y;
        lam_pos(i) = colm;

        % updating the residual
        r = y - A_new*lambda;
        
        % increase increament
        i = i+1;
               
        if(vecnorm(r) <= 0.01)
            break;
        end
    end
    a_x(lam_pos) = lambda;
end


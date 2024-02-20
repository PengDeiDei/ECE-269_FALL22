% Noisy Case
% The sparsity is known
%% N = 20
N = 20;
M = N;
S_max = N;
std_sm = 0.01;
std_lg = 0.1;
success_prob_sm_20 = zeros(S_max,M);
avg_error_sm_20 = zeros(S_max,M);
success_prob_lg_20 = zeros(S_max,M);
avg_error_lg_20 = zeros(S_max,M);

for m = 1:N
    for s = 1:S_max
        success_sm = 0;
        success_lg = 0;

        err_sm = 0;
        err_lg = 0;

        for i = 1:2000
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
            
            % generate noise followed normal distrubution with 0 mean
            n_sm = randn(m,1)*std_sm;
            n_lg = randn(m,1)*std_lg;

            x = zeros(1,N); % generate sparse vector
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';

            y_sm = A*x'+n_sm;
            y_lg = A*x'+n_lg;

            x_approx_sm = OMP(y_sm,A,s);
            x_approx_lg = OMP(y_lg,A,s);

             if(vecnorm(x_approx_sm-x)/vecnorm(x) <= 10^-3)
                 success_sm = success_sm+1;
             else
                 err_sm = err_sm+1;
             end

             if(vecnorm(x_approx_lg-x)/vecnorm(x) <= 10^-3)
                 success_lg = success_lg+1;
             else
                 err_lg = err_lg+1;
             end
        end
        success_prob_sm_20(s,m) = success_sm;
        avg_error_sm_20(s,m) = err_sm;

        success_prob_lg_20(s,m) = success_lg;
        avg_error_lg_20(s,m) = err_lg;
    end
end
%%
figure(1)
imagesc(success_prob_sm_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 20, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(2)
imagesc(success_prob_lg_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 20, std = 0.1']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(3)
imagesc(avg_error_sm_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 20, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(4)
imagesc(avg_error_sm_20/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 20, std = 0.1']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

%% N = 50
N = 50;
M = N;
S_max = N;
std_sm = 0.01;
std_lg = 0.1;
success_prob_sm_50 = zeros(S_max,M);
avg_error_sm_50 = zeros(S_max,M);
success_prob_lg_50 = zeros(S_max,M);
avg_error_lg_50 = zeros(S_max,M);

for m = 1:N
    for s = 1:S_max
        success_sm = 0;
        success_lg = 0;

        err_sm = 0;
        err_lg = 0;

        for i = 1:2000
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
            
            % generate noise followed normal distrubution with 0 mean
            n_sm = randn(m,1)*std_sm;
            n_lg = randn(m,1)*std_lg;

            x = zeros(1,N); % generate sparse vector
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';

            y_sm = A*x'+n_sm;
            y_lg = A*x'+n_lg;

            x_approx_sm = OMP(y_sm,A,s);
            x_approx_lg = OMP(y_lg,A,s);

             if(vecnorm(x_approx_sm-x)/vecnorm(x) <= 10^-3)
                 success_sm = success_sm+1;
             else
                 err_sm = err_sm+1;
             end

             if(vecnorm(x_approx_lg-x)/vecnorm(x) <= 10^-3)
                 success_lg = success_lg+1;
             else
                 err_lg = err_lg+1;
             end
        end
        success_prob_sm_50(s,m) = success_sm;
        avg_error_sm_50(s,m) = err_sm;

        success_prob_lg_50(s,m) = success_lg;
        avg_error_lg_50(s,m) = err_lg;
    end
end
%%
figure(5)
imagesc(success_prob_sm_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 20, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar
figure(6)
imagesc(success_prob_lg_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 50, std = 0.1']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(7)
imagesc(avg_error_sm_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 50, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(8)
imagesc(avg_error_sm_50/2000);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 50, std = 0.1']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

%% N = 100
N = 100;
M = N;
S_max = N;
std_sm = 0.01;
std_lg = 0.1;
success_prob_sm_100 = zeros(S_max,M);
avg_error_sm_100 = zeros(S_max,M);
success_prob_lg_100 = zeros(S_max,M);
avg_error_lg_100 = zeros(S_max,M);

for m = 1:N
    for s = 1:S_max
        success_sm = 0;
        success_lg = 0;

        err_sm = 0;
        err_lg = 0;

        for i = 1:200
            % generate matrix A
            A_mat = randn(m,N);
            A = A_mat./vecnorm(A_mat); % normalization
            
            % generate noise followed normal distrubution with 0 mean
            n_sm = randn(m,1)*std_sm;
            n_lg = randn(m,1)*std_lg;

            x = zeros(1,N); % generate sparse vector
            % asign s random entries within the interval [1,10] to sparse vector
            sign = randi([0 1],s,1);
            x(randperm(N,s)) = randi([1,10],1,s).*((-1).^sign)';

            y_sm = A*x'+n_sm;
            y_lg = A*x'+n_lg;

            x_approx_sm = OMP(y_sm,A,s);
            x_approx_lg = OMP(y_lg,A,s);

             if(vecnorm(x_approx_sm-x)/vecnorm(x) <= 10^-3)
                 success_sm = success_sm+1;
             else
                 err_sm = err_sm+1;
             end

             if(vecnorm(x_approx_lg-x)/vecnorm(x) <= 10^-3)
                 success_lg = success_lg+1;
             else
                 err_lg = err_lg+1;
             end
        end
        success_prob_sm_100(s,m) = success_sm;
        avg_error_sm_100(s,m) = err_sm;

        success_prob_lg_100(s,m) = success_lg;
        avg_error_lg_100(s,m) = err_lg;
    end
end
%%
figure(9)
imagesc(success_prob_sm_100/200);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 100, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(10)
imagesc(success_prob_lg_100/200);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Probability ' ...
    'of Exact Support Recovery, N = 100, std = 0.1']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(11)
imagesc(avg_error_sm_100/200);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 100, std = 0.01']);
xlabel('S_{max}');
ylabel('Measuremnt M');
colorbar

figure(12)
imagesc(avg_error_sm_100/200);
ax = gca;
ax.YDir = 'normal';
title(['The Noiseless Phase Transition of The Average ' ...
    'Normalized Error, N = 100, std = 0.1']);
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
               
        if(vecnorm(r) <= 0.05)
            break;
        end
    end
    a_x(lam_pos) = lambda;
end

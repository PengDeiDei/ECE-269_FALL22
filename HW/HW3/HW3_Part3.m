% Decode Compressed Message
load("Y1 Y2 Y3 and A1 A2 A3.mat");
%% Simply display the compressed images
figure(4)
subplot(2,2,1);
image(y1);
title('The Compressed Image y1');
subplot(2,2,2);
image(y2);
title('The Compressed Image y2');
subplot(2,2,3);
image(y3);
title('The Compressed Image y3');
%% Decode the compressed images using OMP
a_x1 = OMP(y1,A1,size(A1,2));

a_x2 = OMP(y2,A2,size(A2,2));

a_x3 = OMP(y3,A3,size(A3,2));

%% Reshape and display
A_x1 = reshape(a_x1,[90,160]);
A_x2 = reshape(a_x2,[90,160]);
A_x3 = reshape(a_x3,[90,160]);

%%
figure(1)
imagesc(A_x1);
title('The Decoded Image of Input y1');
figure(2)
imagesc(A_x2);
title('The Decoded Image of Input y2');
figure(3)
imagesc(A_x3);
title('The Deocded Image of Input y3');

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


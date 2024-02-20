clc;clear;
% Load images for training
% All neutral expression images are reshaped from (M,N) to (M*N,1)
img = imread("Data for HW4\upload_dataset\2a.jpg");
img_size = size(img);
imgs = zeros([img_size(1)*img_size(2),171]);
idx = 1;
for i = 2: 200
    img_name = "Data for HW4\upload_dataset\"+int2str(i)+"a.jpg";
    if isfile(img_name)
        imgs(:,idx) = reshape(imread(img_name),[img_size(1)*img_size(2),1]);
        idx = idx+1;
    end
end

% Seperate the 1-D images vectors into sets of training and validation
trains = imgs(:,1:100);
valids = imgs(:,101:171);

% Calculate the average of all images in training set
Psi = sum(imgs,2)/size(imgs,2);
% Normalization
A_trains = trains - Psi;
A_valids = valids - Psi;
% Calculate C' instad of C to reduce computation, C' has the same
% eigenvalue as C has
C_p = A_trains'*A_trains;

% Implement SVD on C'
[U,S,V] = svd(C_p);

% *************************** Part 1 *************************************
% Find the singular value of data in training set
singulars = diag(S);


% Calculate the eigenfaces by U_i = A*v_i
eig_faces = A_trains*V;
% Find the 10 most representative eigenfaces 
[~,idx] = sort(sum(eig_faces),'descend');
eig_faces_10 = eig_faces(:,idx(1:10));
% Normalize
eig_faces_10 = eig_faces_10./norm(eig_faces_10);

% Plot the singular value of data in training set
figure
plot(1:100,singulars,'-o','LineWidth',1.5,'MarkerSize',3);
title('The Singular Value of the Neutral Expression Training Data');
ylabel('Singular Value');
grid on;

% Show the 10 most representative eigenfaces
figure
for i = 1:10
    subplot(3,4,i)
    imagesc(reshape(eig_faces_10(:,i),img_size));
    colormap(gray(255));
    subtitle("Eigenface: #"+int2str(idx(i)));
end
sgtitle(['The Most 10 Representative Eigenfaces in ' ...
    'Training Set of Neutral Expression']);
%%
% Randomly pick a normalized image from tranining set of neutral expression
test_id = randi(100);
img_test = A_trains(:,test_id);
% Calculate the weights of each of the 10 max eigenfaces for image reconstruction
w_10 = zeros(10,1);
for i = 1:10
    w_10(i) = eig_faces_10(:,i)'*(img_test); 
end
img_re = eig_faces_10*w_10;
err = abs(sum((img_test - img_re).^2)/(img_size(1)*img_size(2)));

% Display the reconstructed image with MSE
figure
subplot(1,2,1)
imagesc(reshape(img_test+Psi,img_size));
colormap(gray(255));
title("Original Train Image #" + test_id);
subplot(1,2,2)
imagesc(reshape(img_re+Psi,img_size));
colormap(gray(255));
title(["Reconstructed Image using "+length(w_10)+" max eigenfaces", ...
    "error = " + err]); 
%%
% ******************************* Part 2 *********************************
% Reconstruct the chosen image in part 1 with 50 PCs
[img_re, w_n, err_n ]= W_N(50,eig_faces,img_size,img_test);

figure
subplot(1,2,1)
imagesc(reshape(img_test+Psi,img_size));
colormap(gray(255));
title("Original Train Image #" + test_id);
subplot(1,2,2)
imagesc(reshape(img_re+Psi,img_size));
colormap(gray(255));
title(["Reconstructed Image using "+length(w_n)+" max eigenfaces", ...
    "error = " + err]); 

% Calculate and plot the MSE vs. number of PCs
MSEs = zeros(100,1);

for i = 1: 100
    [~,~,MSEs(i)] = W_N(i,eig_faces,img_size,img_test);
end

figure
plot(1:100,MSEs,'-o','LineWidth',1.5,'MarkerSize',3);
title("The MSE vs. number of PCs used of Train Image #"+test_id);
ylabel('MSE');
xlabel('Number of PCs')
grid on;
%%
% ****************************** Part 3 **********************************
% Randomly pick an image from validation set
valid_id = randi(71);
img_valid = A_valids(:,valid_id);

MSEs_valid = zeros(100,1);
for i = 1: 100
    [~,~,MSEs_valid(i)] = W_N(i,eig_faces,img_size,img_valid);
end

figure
plot(1:100,MSEs_valid,'-o','LineWidth',1.5,'MarkerSize',3);
title("The MSE vs. number of PCs used of Validation Image #"+valid_id);
ylabel('MSE');
xlabel('Number of PCs')
grid on;

%% From the plot of MSE vs # of PCs, I choose #PCs = 90 to reconstruct the image
[img_valid_re,~,err_val] = W_N(90,eig_faces,img_size,img_valid);
figure
subplot(1,2,1)
imagesc(reshape(img_valid+Psi,img_size));
colormap(gray(255));
title("Original Validation Image #"+valid_id);
subplot(1,2,2)
imagesc(reshape(img_valid_re+Psi,img_size));
colormap(gray(255));
title(["Reconstructed Image using 90 max eigenfaces", ...
    "error = " + err_val]); 
%%
% Function W_N uses N number of PCs to reconstruct the input image by N PCs,
% and return the weights of PCs and the MSE error as well.
function [img_re, w_n, err]= W_N(N,eig_faces,img_size,img)
    [~,idx] = sort(sum(eig_faces),'descend');
    eig_faces_n = eig_faces(:,idx(1:N));
    eig_faces_n = eig_faces_n/norm(eig_faces_n);
    w_n = zeros(N,1);
    for i = 1:N
        w_n(i) = eig_faces_n(:,i)'*(img); 
    end
    img_re = eig_faces_n*w_n;
    err = abs(sum((img - img_re).^2)/(img_size(1)*img_size(2)));    
end

clc;clear;
% Load images for training
% All smiling images are reshaped from (M,N) to (M*N,1)
smile = imread("Data for HW4\upload_dataset\2b.jpg");
s_size = size(smile);
smiles = zeros([s_size(1)*s_size(2),171]);
neu = imread("Data for HW4\upload_dataset\2a.jpg");
n_size = size(neu);
neus = zeros([n_size(1)*n_size(2),171]);

idx = 1;
for i = 2: 200
    smile_name = "Data for HW4\upload_dataset\"+int2str(i)+"b.jpg";
    if isfile(smile_name)
        smiles(:,idx) = reshape(imread(smile_name),[s_size(1)*s_size(2),1]);
        idx = idx+1;
    end
end

idx = 1;
for i = 2: 200
    neu_name = "Data for HW4\upload_dataset\"+int2str(i)+"a.jpg";
    if isfile(neu_name)
        neus(:,idx) = reshape(imread(neu_name),[n_size(1)*n_size(2),1]);
        idx = idx+1;
    end
end

% Smiling set 
% Seperate the 1-D images vectors into sets of training and validation
s_trains = smiles(:,1:100);
s_valids = smiles(:,101:171);
% Calculate the average of all smiling images in training set
Psi_s = sum(smiles,2)/size(smiles,2);
% Normalization
s_trains = s_trains - Psi_s;
s_valids = s_valids - Psi_s;
% Calculate C' instad of C to reduce computation, C' has the same
% eigenvalue as C has
C_p_s = s_trains'*s_trains;
% Implement SVD on C' of smiling set
[~,~,V_s] = svd(C_p_s);
% Calculate the eigenfaces by U_i = A*v_i
eig_faces_s = s_trains*V_s;

% Neutral Set 
n_trains = neus(:,1:100);
n_valids = neus(:,101:171);

Psi_n = sum(neus,2)/size(neus,2);

n_trains = n_trains - Psi_n;
n_valids = n_valids - Psi_n;

C_p_n = n_trains'*n_trains;

[~,~,V_n] = svd(C_p_n);
eig_faces_n = n_trains*V_n;

%% In part 1, I choose to use 10 max eigenfaces to reconstruct the image
% Combine both two validations set of smiling and neutral expression images
tests_all = [s_valids, n_valids];
% Randomly pick 60 images from them
test_ids = randperm(size(tests_all,2),60);
tests = tests_all(:,test_ids);
% Generate the labels of the test set
label = zeros(60,1);
for i = 1:60
    id = test_ids(i);

    if id > 71
        label(i) = 1;
    end

end

preds = ones(60,1)*(-1);
acc_s = 0;
acc_n = 0;

for i = 1:60
    [s_re,~,err_s] = W_N(10,eig_faces_s,s_size,tests(:,i));
    [n_re,~,err_n] = W_N(10,eig_faces_n,n_size,tests(:,i));

    if err_s>err_n
        pred = 0;
    else
        pred = 1;
    end

    if pred == labels(i)
       if pred == 0
           acc_s = acc_s+1;
       else
           acc_n = acc_n+1;
       end
    else
        figure
        subplot(1,3,1)
        imagesc(reshape(tests(:,i),s_size));
        colormap(gray(255));
        title("Missed Original Test Image #"+i);
        subplot(1,3,2)
        imagesc(reshape(s_re,s_size));
        colormap(gray(255));
        title(["Reconstructed Smiling Image #"+i, err_s]);
        subplot(1,3,3)
        imagesc(reshape(n_re,n_size));
        colormap(gray(255));
        title(["Reconstructed Neutral Image #"+i,err_n]);
    end
end
%%
ac_n = acc_n/60;
ac_s = acc_s/60;
miss = 1-(ac_s+ac_n);
disp("The accurate rate of neutral expression prediction = "+ac_n);
disp("The accurate rate of smiling expression prediction = "+ac_s);
disp("The miss rate of prediction = "+miss);
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

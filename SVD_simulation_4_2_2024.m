clear all
clc
close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

%%
% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\toluene_prompt\';
% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Cycloheptatriene_prompt\';
% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_5fs';
% data = load_csv_data(directory_path);
load("data_toluene.mat")

frag_name_list = {'C','C','C','C','C','C','C','H','H','H','H','H','H','H','H'};
for i = 1:15;

    data.frag_param{i}.name = frag_name_list{i};

end

%%

% struct = cat(3,data.CoM_mom{1},data.CoM_mom{2},data.CoM_mom{3},data.CoM_mom{4},data.CoM_mom{5},data.CoM_mom{6},data.CoM_mom{7});
nPoints = length(data.CoM_mom{1}(:,3));
X = [data.CoM_mom{1}(:,3),data.CoM_mom{2}(:,3)]';
Xavg = mean(X,2);
B = X-Xavg*ones(1,nPoints);
[U,S,V] = svd(X/sqrt(nPoints),'econ'); % PCA via SVD
theta = (0:0.1:1)*2*pi;
Xstd = U*S*[cos(theta);sin(theta)];
setup
figure
det_image_plot(X(1,:),X(2,:),-200:5:200,-200:5:200)
hold on
plot(Xavg(1)+Xstd(1,:),Xavg(2) + Xstd(2,:),'b--',"LineWidth",5)

% u1 = U(:,1);
% u2 = U(:,2);
% sigma1 = S(1,1);
% sigma2= S(2,2);

% figure
% det_image_plot(X(:,1),X(:,2),-200:5:200,-200:5:200)
% hold on
% plot(sigma1*u1)
% plot(sigma2*u2)


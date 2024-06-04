clear all
clc
close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\toluene_prompt\';
% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Cycloheptatriene_prompt\';
% directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_5fs';
% data = load_csv_data(directory_path);
load("data_toluene.mat")
data.XYT = data.CoM_mom;
frag_name_list = {'C','C','C','C','C','C','C','H','H','H','H','H','H','H','H'};
for i = 1:15;

    data.frag_param{i}.name = frag_name_list{i};

end

for j = 1:length(data.CoM_mom{i})
    R_j = RandomRotationMatrix;
    for i = 1:15;
        data.CoM_mom{i}(j,:) = R_j*data.CoM_mom{i}(j,:)';
    end
end

data = rand_perm(data,[1 2 3 4 5 6 7]);
% data = rand_perm(data,[8 9 10 11 12 13 14 15]);

data = false_coin(data,0.0);


figure; 
for i =1:3;
    histogram(sqrt(sum(data.CoM_mom{i}.^2,2)),0:1:300,"DisplayStyle","stairs","LineWidth",2)
    hold on
end
xlabel('|p_{i}|^{2} a.u')

plot_tot(data)

data = mol_frame_calc(data,1,2,[3],1000000000);

figure
[tx,ty,tz] = projection_maps(data.Mol_mom{3}(:,1),data.Mol_mom{3}(:,2),data.Mol_mom{3}(:,3),...
-3:0.05:3,false);


%%
tic

shots = 9999;


cmap = turbo;
clr = cmap(round(linspace(1,length(cmap),21)),:);
p_ij_norm_all = [];
p_ij_ang_all = [];
q=1;

p_ij_norm_edges = 200:6:600;
p_ij_norm_bins = md_pts(p_ij_norm_edges);
p_ij_ang_edges = 0:4:180;
p_ij_ang_bins = md_pts(p_ij_ang_edges);
p_ij_dataset = zeros(shots,length(p_ij_norm_bins),length(p_ij_ang_bins));

% shots = 9999;
for shot = 1:shots;    
    p_ij_norm_shot = [];
    p_ij_ang_shot = [];
    for i = 1:3
        for j=1:3
            if i<j;
                p_ij_norm = sqrt(sum((data.CoM_mom{i}(shot,:)-data.CoM_mom{j}(shot,:)).^2,2));
                p_ij_ang = vec_angle(data.CoM_mom{i}(shot,:),data.CoM_mom{j}(shot,:));
                p_ij_norm_shot = [p_ij_norm_shot; p_ij_norm];
                p_ij_ang_shot = [p_ij_ang_shot ; p_ij_ang];


%                 p_ij_norm_all = [p_ij_norm_all; p_ij_norm];
%                 p_ij_ang_all = [p_ij_ang_all; p_ij_ang];
    %             scatter(p_ij_norm,p_ij_ang,'MarkerFaceColor',clr(q,:),"MarkerEdgeColor","none","DisplayName",[num2str(i) ' - ' num2str(j)]);
    %             subplot(3,7,q)
    %             det_image_plot(p_ij_norm,p_ij_ang,100:1:800,0:1:180);
    %             q=q+1;
    %             title([num2str(i) ' - ' num2str(j)])


            end
        end
    end

    Nxy = histcounts2(p_ij_norm_shot,p_ij_ang_shot,p_ij_norm_edges,p_ij_ang_edges);
    p_ij_dataset(shot,:,:) = Nxy;


end
% legend
% fontsize(gcf, 10, 'points')

toc
%%

figure

% det_image_plot(p_ij_norm_all,p_ij_ang_all,200:1:600,0:1:180);

p_ij_dataset_shot_int = reshape(sum(p_ij_dataset),size(p_ij_dataset,2),size(p_ij_dataset,3));
pcolor(p_ij_norm_bins,p_ij_ang_bins,p_ij_dataset_shot_int');
shading flat


xlabel("||p_{i} - p_{j}||")
ylabel("\theta_{ij}")
colorbar

%%

tic
p_ij_dataset_reshape = reshape(p_ij_dataset,prod(size(p_ij_dataset,[2,3])),size(p_ij_dataset,1));
[U,S,V] = svd(p_ij_dataset_reshape);
% toc

%%

cov_dataset = cov(p_ij_dataset_reshape');


%%
setup
tic
S_red = zeros(size(S,1),size(S,2));
S_red(1,1) = S(1,1);
p_ij_dataset_red = reshape(U*S_red*V',prod(size(p_ij_dataset,1)),size(p_ij_dataset,2),size(p_ij_dataset,3));
p_ij_dataset_red_int = reshape(sum(p_ij_dataset_red),size(p_ij_dataset,2),size(p_ij_dataset,3));
toc

figure
subplot(2,2,1)
semilogy(diag(S),'*-')
hold on
semilogy(diag(S_red),'*-')
xlabel("Index")
ylabel("SVD")

subplot(2,2,2)
pcolor(p_ij_norm_bins,p_ij_ang_bins,p_ij_dataset_shot_int');
shading flat

xlabel("||p_{i} - p_{j}||")
ylabel("\theta_{ij}")
title('Raw')
subplot(2,2,3)
pcolor(p_ij_norm_bins,p_ij_ang_bins,p_ij_dataset_red_int');
shading flat

xlabel("||p_{i} - p_{j}||")
ylabel("\theta_{ij}")
title('Recon')

%%

file = {"no_false_coin.mat","50_false_coin.mat","100_false_coin.mat"};
name = {0,50,100};
figure
for i = 1:length(file);
    load(file{i})
    semilogy(S_diag,"*-","DisplayName",num2str(name{i}))
    hold on
end
legend

%%

function plot_tot(data);
    
    tot = zeros(length(data.CoM_mom{1}),3);
    
    for i = 1:7;
        tot = tot + data.CoM_mom{i};
    end
    
    bins = -1000:50:1000;

    figure
    for j =1:3;
    histogram(tot(:,j),bins,"DisplayStyle","stairs","LineWidth",3);
    hold on

    end
    xlabel('Momentum Sum / a.u')

end



function output = false_coin(data,percentageToSwap)

    
    output = data;
    
    % Define the percentage of rows to swap (10%)
    
    % Iterate through each entry in the cell array
    for i = 1:7
        % Get the matrix from the current cell
        matrix = output.CoM_mom{i};
        
        % Determine the number of rows to swap (10% of total rows)
        numRows = size(matrix, 1);
        numToSwap = round(percentageToSwap * numRows);
        
        % Randomly select rows to swap
        rowsToSwap = randperm(numRows, numToSwap);
        
        % Randomly select rows to swap with
        rowsToSwapWith = randperm(numRows, numToSwap);
        
        % Swap selected rows
        for j = 1:numToSwap
            % Store the row to be replaced
            tempRow = matrix(rowsToSwap(j), :);
            % Perform the swap
            matrix(rowsToSwap(j), :) = matrix(rowsToSwapWith(j), :);
            matrix(rowsToSwapWith(j), :) = tempRow;
        end
        
        % Update the matrix in the cell array
        output.CoM_mom{i} = matrix;
    end

    rand_id = randperm(length(data.CoM_mom{1}(:,1)));

    for i = 1:15
        output.CoM_mom{i} = output.CoM_mom{i}(rand_id,:);
    end


end
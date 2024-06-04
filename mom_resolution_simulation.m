clear all
clc
close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

%%
directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_prompt\';
% % directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_5fs';
% data = load_csv_data(directory_path);
load("data_toluene.mat")


N = 100000;
%% Get single geometry

frag_name_list = {'C','C','C','C','C','C','C','H','H','H','H','H','H','H','H'};
for i = 1:15;
    data.frag_param{i}.name = frag_name_list{i};
    data.CoM_mom_start{i} = repmat(data.CoM_mom{i}(1,:),N,1);
end

for n = 1:N;

% [Q,R] = qr(randn(3));
[Q,~] = qr(randn(3));Q(:,1)=Q(:,1)*(2*(rand>0.5)-1);Q(:,2)=det(Q)*Q(:,2);

for i = 1:15
    data.CoM_mom_start{i}(n,:) = Q*data.CoM_mom_start{i}(n,:)';
end

end

%% Genearte normal distribution of the momentum for each fragment

p_res = 0.1;
for i = 1:15
    
    data.CoM_mom{i} = data.CoM_mom_start{i} + normrnd(0,p_res*abs(data.CoM_mom_start{i}));

end


data_mat = cell2mat(data.CoM_mom);

mom_sum_bins = -200:5:200;
figure
histogram(sum(data_mat(:,1:3:45),2),mom_sum_bins,"DisplayName",'x')
hold on
histogram(sum(data_mat(:,2:3:45),2),mom_sum_bins,"DisplayName",'y')
histogram(sum(data_mat(:,3:3:45),2),mom_sum_bins,"DisplayName",'z')
title(['Percent Resolution:  ' num2str(100*p_res) '%'])
legend
%%

card_dir = {'x','y','z'};
q=1; qn = card_dir{q};
r=2; rn = card_dir{r};
figure
p_bin = -1.25:0.025:1.25;

plot_indices = reshape(1:7*7,7,7)';
% t = tiledlayout(7,7);

for i = 1:7;
    for j = 1:7;
        
        if i ~= j;

            plt_index = [1,2,3,4,5,6,7];
            plt_index([i j]) = [];
            data = mol_frame_calc(data,i,j,plt_index,1000000000);
            
            px = []; py = []; pz = [];
            for k = plt_index;
            
                
                px= [px; data.Mol_mom{k}(:,1)];
                py= [py; data.Mol_mom{k}(:,2)];
                pz= [pz; data.Mol_mom{k}(:,3)];
            
            
            end
            p_dir{1} = px; p_dir{2} = py; p_dir{3} = pz;
            
            subplot(7,7,plot_indices(i,j))
%             histogram(px,p_bin,"DisplayStyle","stairs")
%             hold on
%             histogram(py,p_bin,"DisplayStyle","stairs")
%             histogram(pz,p_bin,"DisplayStyle","stairs","LineWidth",2.5,"EdgeColor",[0.4941    0.1843    0.5569])
            subplot(7,7,plot_indices(i,j))
            
            det_image_plot(p_dir{q},p_dir{r},p_bin,p_bin);
            ax = gca;
            pos = ax.Position;
            axis square
%             axis off
            ax.Position= [pos(1) pos(2) 0.1 0.1];
%             ax.FontSize = 1; 

        else

            ax = gca;
%             set(ax,'XTick',[], 'YTick', []);        

        end
    end
end

sgtitle([qn ' vs ' rn])

% save_to_clipboard


%%

plot_indices = reshape(1:3*3,3,3)';


px = []; py = []; pz = [];

for i = 1:7;
    for j = 1:7;
        
        if i ~= j;
            
%             plt_index = [8,9,10,11,12,13,14,15];
            plt_index = [1,2,3,4,5,6,7];
            plt_index([i j]) = [];
            data = mol_frame_calc(data,i,j,plt_index,1000000000);
            
            for k = plt_index;
            
                
                px= [px; data.Mol_mom{k}(:,1)];
                py= [py; data.Mol_mom{k}(:,2)];
                pz= [pz; data.Mol_mom{k}(:,3)];
            
            
            end
            p_dir{1} = px; p_dir{2} = py; p_dir{3} = pz;
            
            
%             det_image_plot(p_dir{q},p_dir{r},p_bin,p_bin);
%             subplot(3,3,plot_indices(i,j))
%             projection_maps(px,py,pz,...
%                 -3:.01:3);

        end
    end
end



%%
figure
projection_maps(px,py,pz,...
-3:.0005*5:3);


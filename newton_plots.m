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
%                 px= [data.Mol_mom{k}(:,1)];
%                 py= [data.Mol_mom{k}(:,2)];
%                 pz= [data.Mol_mom{k}(:,3)];            
            
            end
            p_dir{1} = px; p_dir{2} = py; p_dir{3} = pz;


        end
    end


end



%%


figure
[tx,ty,tz] = projection_maps(px,py,pz,...
-3:.025/4:3,false);
% title(['i= ' num2str(i) '  j= ' num2str(j)])


%%
p_bins = -1.5:0.025/2:1.5;

jet_0 = jet;
color_list = jet_0(round(linspace(1,length(jet_0),5)),:);
setup

for i = 1:7;
    for j = 1:7;


        
        if i ~= j;


figure
subplot(2,2,1)
hold on;
det_image_plot(px,py,p_bins,p_bins)
xlabel('x'); ylabel('y');axis square

subplot(2,2,2)
hold on;
det_image_plot(py,pz,p_bins,p_bins)
xlabel('y'); ylabel('z');axis square

subplot(2,2,3)
hold on;
det_image_plot(pz,px,p_bins,p_bins)
xlabel('z'); ylabel('x');axis square

            
%             plt_index = [8,9,10,11,12,13,14,15];
            plt_index = [1,2,3,4,5,6,7];
            plt_index([i j]) = [];
            data = mol_frame_calc(data,i,j,plt_index,1000000000);
            
            q= 1;
            for k = plt_index;
            
                
                px_k= [data.Mol_mom{k}(:,1)];
                py_k= [data.Mol_mom{k}(:,2)];
                pz_k= [data.Mol_mom{k}(:,3)];            

                subplot(2,2,1)
                pt_xy = errorbar(mean(px_k),mean(py_k),std(py_k),std(py_k),std(px_k),std(px_k),...
                    'o',"LineWidth",1.5,"MarkerFaceColor",color_list(q,:),"Color",color_list(q,:),"MarkerEdgeColor",'w');
                subplot(2,2,2)
                pt_yz = errorbar(mean(py_k),mean(pz_k),std(pz_k),std(pz_k),std(py_k),std(py_k),...
                    'o',"LineWidth",1.5,"MarkerFaceColor",color_list(q,:),"Color",color_list(q,:),"MarkerEdgeColor",'w');
                subplot(2,2,3)
                pt_zx = errorbar(mean(pz_k),mean(px_k),std(px_k),std(px_k),std(pz_k),std(pz_k),...
                    'o',"LineWidth",1.5,"MarkerFaceColor",color_list(q,:),"Color",color_list(q,:),"MarkerEdgeColor",'w');

                q=q+1;

            end



        end
    end


end


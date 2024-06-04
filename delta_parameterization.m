clear all
clc
close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_1fs')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_2fs')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_5fs')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_10fs')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_100fs')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

%%
directory_path = 'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_1fs';
data = load_csv_data(directory_path);
% load("data_toluene.mat")

frag_name_list = {'C','C','C','C','C','C','C','H','H','H','H','H','H','H','H'};
for i = 1:15;

    data.frag_param{i}.name = frag_name_list{i};

end

%%
main_path = {'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_1fs\',...
    'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_5fs\',...
    'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_charge_build_v3\delta_10fs\',...
    'C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8 Geometry Files\Toluene_prompt\'};

figure

for l = 1:length(main_path)
        
    directory_path = main_path{l};
%     dir_list_l = dir_list{l};
%     directory_path = [dir_list_l];
    if l == 4;
        load('data_toluene.mat')
    else
        data = load_csv_data(directory_path);
    end
    
    for i = [1 2 3 4 5 6 7]
    subplot(3,3,i)
        ke_i = dot(data.CoM_mom{i},data.CoM_mom{i},2)/(2*12*1822.88)*27;
        ke_i = ke_i(1:500);
    
        histogram(ke_i,0:0.5:60);
        hold on
        
        xlabel('KE')
        title(['C: ' num2str(i)])
    end

end
    legend('1 fs', '5 fs','10 fs', 'prompt')




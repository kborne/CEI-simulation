clear all
clc
close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

load("data_toluene.mat")

frag_name_list = {'C','C','C','C','C','C','C','H','H','H','H','H','H','H','H'};
for i = 1:15;

    data.frag_param{i}.name = frag_name_list{i};

end

%%


card_dir = {'x','y','z'};
q=1; qn = card_dir{q};
r=2; rn = card_dir{r};

p_bin = -1.25:0.025:1.25;
plot_indices = reshape(1:7*7,7,7)';
i=1; j= 3;


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
  
%%

figure
tst = scatter(px,py);
sav


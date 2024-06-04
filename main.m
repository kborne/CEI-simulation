clear all
clc
% close all

addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448\functions')
addpath('C:\Users\kborn\Analysis\SQS-EuXFEL\C7H8-2448')

%%


% load("data_toluene.mat")
% load('data_cyclo.mat')
% load("data_hep.mat")
load('data_hep_GG_trans.mat')
data_GG_trans = data;
clear data;
load('data_hep_AG.mat')
data_AG = data;
clear data;
for i = 1:15;
    data_i_GG = data_GG_trans.CoM_mom{i};
    data_i_AG = data_AG.CoM_mom{i};
    data.CoM_mom{i} = [data_i_GG; data_i_AG];
end

cmap = jet;
c_list = cmap(round(linspace(1,length(cmap),7)),:);
%%

figure
for i = 1:7;
    pi = data.CoM_mom{i};
    histogram(vecnorm(pi,2,2).^2/(2*12*1822),0:0.01:5,"DisplayName",num2str(i),"DisplayStyle","stairs","LineWidth",3);
    hold on
end
legend
xlabel("Ion Kinetic Energy / a.u.")
ylabel("Counts")
%%




%%
close all
card_dir = {'x','y','z'};
q=2; qn = card_dir{q};
r=3; rn = card_dir{r};

p_bin = -1.25:0.01:1.25;
t = tiledlayout(7,7);

for i = [1,2,3,4,5,6,7];
    for j = [1,2,3,4,5,6,7];
        
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
            

                nexttile
                det_image_plot(p_dir{q},p_dir{r},p_bin,p_bin);
                xline(-1:0.5:1,':')
                yline(-1:0.5:1,':')
                ax = gca;
            axis square    
            set(ax,'XTick',[], 'YTick', []);

        else
            nexttile
            ax = gca;
            axis square
            axis off

        end
    end

end

t.TileSpacing = 'none';


%%


close all
card_dir = {'x','y','z'};
q=2; qn = card_dir{q};
r=3; rn = card_dir{r};

p_bin = -1.25:0.01:1.25;
t = tiledlayout(7,7);

cmap_slct = jet;
clr_list = cmap_slct(round(linspace(1,length(cmap_slct),7)),:);

for i = [1,2,3,4,5,6,7];
    for j = [1,2,3,4,5,6,7];
        
        if i ~= j;

            plt_index = [1,2,3,4,5,6,7];
            plt_index([i j]) = [];
            data = mol_frame_calc(data,i,j,plt_index,1000000000);
            
            nexttile
                
            for k = plt_index;
                        
                px= data.Mol_mom{k}(:,1);
                py= data.Mol_mom{k}(:,2);
                pz= data.Mol_mom{k}(:,3);
                scatter(px,py,'o',"MarkerEdgeColor","none","MarkerFaceColor",clr_list(k,:),"MarkerFaceAlpha",0.25);
                hold on
                text(mean(px),mean(py),num2str(k))

            end
            

                xline(-1:0.5:1,':')
                yline(-1:0.5:1,':')
                ax = gca;
            axis square    
            set(ax,'XTick',[], 'YTick', []);

        else
            nexttile
            ax = gca;
            axis square
            axis off

        end
    end

end

t.TileSpacing = 'none';

%%
plot_indices = reshape(1:3*3,3,3)';


px = []; py = []; pz = [];
p_indx =[];
cos_theta = [];

i_list=1:7;
j_list = 1:7;
for ii =1:length(i_list);
    i = i_list(ii);
    for jj = 1:length(j_list);
        j = j_list(jj);
        
        if i ~= j;
            
            plt_index = [1,2,3,4,5,6,7];
            plt_index([i,j]) = [];
            data = mol_frame_calc(data,i,j,plt_index,1000000000);
            cos_theta_ij= dot(data.CoM_mom{i},data.CoM_mom{j},2)./(vecnorm(data.CoM_mom{i},2,2).*vecnorm(data.CoM_mom{j},2,2));

            for k = plt_index;
            

                px= [px; data.Mol_mom{k}(:,1)];
                py= [py; data.Mol_mom{k}(:,2)];
                pz= [pz; data.Mol_mom{k}(:,3)];    
                p_indx = [p_indx; k*ones(length(data.Mol_mom{k}(:,3)),1)];

                cos_theta = [cos_theta; cos_theta_ij];

            
            end


        end
    end


end

%%
%
setup
figure
[tx,ty,tz] = projection_maps(px,py,pz,...
-2:0.01:2,false);
xlabel("P_{x} / au")
ylabel("P_{y} / au")
zlabel("P_{z} / au")
cb = colorbar('northoutside','Position',[0.2892    0.7648    0.1500    0.0175],"FontSize",10);

% ylabel(cb,'Counts','FontSize',12)
%%
title(['i = ' num2str(i_list) ' :  j =  ' num2str(j_list)])

% spherical_coordinates([px,py,pz],0.0:0.001:2,3.14159/100,3.14159/100);
% suptitle(['i = ' num2str(i_list) ' :  j =  ' num2str(j_list)])

quantiles = [-1 -0.83  -0.3  0.1 0.65 1];
colours = get(gca, 'colororder');
colours = colours(2:end,:);
colours(3,:) = [1,1,1];


[aximuth,elevation,r] = cart2sph(px,py,pz);

figure; det_image_plot(aximuth,elevation,-3.14159:3.14159/60:3.14159,-3.14159/2:3.14159/60:3.14159/2); 
ax =gca;
    ax.FontSize = 16; 
set(gca,'XTick',-pi:pi/4:pi); 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
set(gca,'YTick',-pi/2:pi/4:pi/2); 
set(gca,'YTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
xlabel('\phi / rad.'); ylabel('\theta / rad.')
cb = colorbar('northoutside',"FontSize",10,'Position',[0.1471+0.55 0.7765+0.1 0.25 0.0434]);
% title("No Gate");

%%

figure; histogram(cos_theta,-1:0.01:1,"DisplayStyle","stairs","EdgeColor",'k'); xlabel('Cos(\theta_{ij})'); ax_hist = gca;hold on;
colours= [[ 1  0  0];
 [ 0.4667    0.6745    0.1882]; 
  [1 1 1];
  [0.4941    0.1843    0.5569];
  [1 1 1]];

figure; 
for i = 1:length(quantiles)-1


    gate = (quantiles(i)<cos_theta) & (cos_theta<quantiles(i+1)); 
    sum(gate/length(gate))

%     histogram(ax_hist,cos_theta(gate),-1:0.01:1,"EdgeColor","none","FaceColor",colours(i,:));
%     figure; projection_maps(px(gate),py(gate),pz(gate),-3:0.01:3,false); title([num2str(quantiles(i:i+1))]);
%     cb = colorbar('northoutside',"FontSize",10,'Position',[0.306311309523807,0.750476190476191,0.125117261904764,0.0182420754911]);
%     title(['cos(\theta_{ij}) = \{' num2str(quantiles(i)) ',' num2str(quantiles(i+1))  '\}']);
%     [aximuth,elevation,r] = cart2sph(px(gate),py(gate),pz(gate));
    
    figure; det_image_plot(aximuth(gate),elevation(gate),-3.14159:3.14159/60:3.14159,-3.14159/2:3.14159/60:3.14159/2);
ax =gca;
    ax.FontSize = 16; 
set(gca,'XTick',-pi:pi/4:pi); 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
set(gca,'YTick',-pi/2:pi/4:pi/2); 
set(gca,'YTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
xlabel('\phi / rad.'); ylabel('\theta / rad.')

cb = colorbar('northoutside',"FontSize",10,'Position',[0.1471+0.55 0.7765+0.1 0.25 0.0434]);
title([num2str(quantiles(i:i+1))]);

end

%%

[aximuth,elevation,r] = cart2sph(px,py,pz);
for i = 1:length(quantiles)-1


gate_cos = (quantiles(i)<cos_theta) & (cos_theta<quantiles(i+1)); 
figure
for k = 3:6
subplot(2,2,k-2);

p_indx_slct = k==p_indx;
gate = p_indx_slct & gate_cos;
det_image_plot(aximuth(gate ),elevation(gate ),-3.14159:3.14159/60:3.14159,-3.14159/2:3.14159/60:3.14159/2); title(['ion = ' num2str(k)]);
xlabel('azimuth'); ylabel('elevation')

end
suptitle([num2str(quantiles(i:i+1))])
end

%%
figure; 
for i = 1:length(quantiles)-1


    gate = (quantiles(i)<cos_theta) & (cos_theta<quantiles(i+1));     
[aximuth,elevation,r] = cart2sph(px(gate),py(gate),pz(gate));  
p_indx_gate = p_indx(gate);

subplot(2,3,i)

title([num2str(quantiles(i:i+1))]);
for k = 3:6;

    p_indx_slct = k==p_indx_gate;
histogram(aximuth(p_indx_slct),-3.14159:3.14159/60:3.1415,"DisplayStyle","stairs","LineWidth",3,"DisplayName",num2str(k),"EdgeColor",clr_list(k,:)); hold on;
end

end
legend

%%
figure
sfh1 = subplot(1,2,1);
sfh1.Position = [0.1347    0.2160    0.3419    0.5907];

plot_indices = reshape(1:7*7,7,7)';
q=1;
cos_theta_tot = [];
cos_theta_avg = zeros(7,7);
redblue_clr = redblue;
for i = 1:7;
    for j = 1:7;
        
        if i~=j;
            p_i = data.CoM_mom{i};
            p_j = data.CoM_mom{j};
            cos_theta_ij= dot(p_i,p_j,2)./(vecnorm(p_i,2,2).*vecnorm(p_j,2,2));
            cos_theta_tot = [cos_theta_tot; cos_theta_ij];
            [~,pt] = min(abs(mean(cos_theta_ij)-linspace(-1,1,length(redblue_clr))))
            histogram(cos_theta_ij,-1:0.01:1,"FaceColor",redblue_clr(pt,:))
            hold on
%             ylim([0,10000])
            cos_theta_avg(i,j)=mean(cos_theta_ij);
        end
        q=q+1
        

    end
end
xlabel('Cos(\theta_{ij})')
ylabel('Counts')
text(-1.5,7000,'E)','FontSize',15)
ax = gca;
ax.FontSize = 16; 

subplot(1,2,2)
C = cos_theta_avg;
image(C,'CDataMapping','scaled');
hold on
clim([-1 1]);
colormap(redblue);
cb = colorbar(); 
ylabel(cb,'Cos(\theta_{ij})','FontSize',16,'Rotation',270)
xlabel('Carbon Number')
ylabel('Carbon Number')


[rows, cols] = size(C);
for i = 1:rows
    for j = 1:cols
        if j~=i
            text(j, i, [num2str(round(C(i, j),2))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'black',"FontSize",15);
            [~,indx] = min(abs(0.5*(quantiles(1:end-1)+quantiles(2:end))-C(i,j)));
        end
%         if j~=i;
%             plot(j,i+0.25,'o','MarkerSize',10','MarkerFaceColor',colours(indx,:),"MarkerEdgeColor",'k')
%         end

        if C(i, j) >= -0.83 && C(i, j) <= -0.3;
            plot(j,i,'sq','MarkerSize',50,'MarkerFaceColor',"none","MarkerEdgeColor",[0.4660 0.6740 0.1880],"LineWidth",7)
        elseif C(i, j) >= 0.1&& C(i, j) <= 0.65;
            plot(j,i,'sq','MarkerSize',50,'MarkerFaceColor',"none","MarkerEdgeColor",[0.4940 0.1840 0.5560],"LineWidth",7)
         elseif C(i, j) >= -1&& C(i, j) <= -0.83;
            plot(j,i,'sq','MarkerSize',50,'MarkerFaceColor',"none","MarkerEdgeColor",[ 1  0  0],"LineWidth",7)
        end
    end
end

axis square
text(-0.5,1,'F)','FontSize',15)
ax = gca;
ax.FontSize = 16; 


%%
pts = 1:7;
plot_slct = 'projection';
for i = 1:7
    for j = 1:7
        if i~=j;
            pts_plt = pts; pts_plt([i,j]) = [];     calc_and_plot_newton(data,i,j,pts_plt,plot_slct)
            saveas(gcf,['newton_maps/cyclo/projection_' 'i=' num2str(i) '_j=' num2str(j)],'png')
        end
    end
end
%%
function [bins] = make_bins(data)
    bin_sz = (max(data)-min(data))/sqrt(length(data));
    bins = min(data):bin_sz:max(data);
end

function [] = calc_and_plot_newton(data_i,x_ref,y_ref,plt_ref,map);

    data_i = mol_frame_calc(data_i,x_ref,y_ref,1:7,10000000);

    px = []; py = []; pz = [];
    for k =1:length(plt_ref);
        px = [px; data_i.Mol_mom{plt_ref(k)}(:,1)];
        py = [py; data_i.Mol_mom{plt_ref(k)}(:,2)];
        pz = [pz; data_i.Mol_mom{plt_ref(k)}(:,3)];
    
    end
    figure
    if strcmp(map,'projection')
        projection_maps(px,py,pz,...
            -2:0.01:2,false)
    elseif strcmp(map,'scatter')
        projection_maps_scatter(data_i.Mol_mom,plt_ref)
    elseif strcmp(map,'angle')
        [aximuth,elevation,r] = cart2sph(px,py,pz);  
        det_image_plot(aximuth,elevation,make_bins(aximuth),make_bins(elevation)); xlabel('azimuth'); ylabel('elevation')
        xlim([-3.14159 3.14159]);
        ylim([-3.14159/2 3.14159/2]);
    elseif strcmp(map,'angle scatter')
        cmap = jet;
        c_list = cmap(round(linspace(1,length(cmap),7)),:);

        for i =1:length(plt_ref);
            px_i = data_i.Mol_mom{plt_ref(i)}(:,1);
            py_i = data_i.Mol_mom{plt_ref(i)}(:,2);
            pz_i = data_i.Mol_mom{plt_ref(i)}(:,3);
            
            [aximuth,elevation,r] = cart2sph(px_i,py_i,pz_i);  
            scatter(aximuth,elevation,'o',"MarkerEdgeColor","none","MarkerFaceColor",c_list(plt_ref(i),:),"MarkerFaceAlpha",0.01); hold on;
            text(mean(aximuth),mean(elevation),num2str(plt_ref(i)))
            xlabel('azimuth'); ylabel('elevation')
        xlim([-3.14159 3.14159]);
        ylim([-3.14159/2 3.14159/2]);    
        
        end
    else
        disp(['no plot called  :  '   map])
    end
    suptitle(['i = ' num2str(x_ref) ' j = ' num2str(y_ref) '  z = ' num2str(plt_ref)]);
end
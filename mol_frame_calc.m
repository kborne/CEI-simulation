function [data_out] = mol_frame_calc(data_in,i_index,j_index,index_plot,pz_gate)

p_i = data_in.CoM_mom{i_index}; 
p_j = data_in.CoM_mom{j_index};

ux = p_i./sqrt(dot(p_i,p_i,2));


uynn = p_j-dot(p_j,ux,2).*ux;  uy = uynn./sqrt(dot(uynn,uynn,2));
uz = cross(uy,ux);


pi_u = [dot(p_i,ux,2),dot(p_i,uy,2),dot(p_i,uz,2)]; n_fac = sqrt(dot(pi_u,pi_u,2));


for k = 1:length(data_in.CoM_mom)  %
    
    pk = data_in.CoM_mom{k};
    pk_u = [dot(pk,ux,2),dot(pk,uy,2),dot(pk,uz,2)]./n_fac;
%     pk_u = [dot(pk,ux,2),dot(pk,uy,2),dot(pk,uz,2)];
    data_in.Mol_mom{k} = pk_u;

    
end


p_xyz_mol  = [];
% gate = sqrt(data.Mol_mom{1}(:,3).^2)<pz_gate | sqrt(data.Mol_mom{2}(:,3).^2)<pz_gate  ;
for i = index_plot


        p_xyz_mol = [p_xyz_mol ; data_in.Mol_mom{i}];
end


data_out = data_in;

end
function [MUA_vec, LFP_vec, VoltVec] = getMUA_soma_neuron(I_mat, I_mat_PS, len_Ts, Nx, Ny, Nz, dx, dy, dz, i, j, y, gamma, sf)

Ve_vec = zeros(len_Ts, 1);

LFP_radius = 0.5;

for t_ind = 1:len_Ts
    
    sum_eP = 0;
    for i_p = 1:Nx
        for j_p = 1:Nz
            
%             if (abs(i_p - i)*dx > LFP_radius) || (abs(j_p - j)*dz > LFP_radius)
%                 continue
%             end
            
            r_s = sqrt(((i_p - i)*dx).^2 + ((j_p - j)*dz).^2 + ((0-y)*dy).^2);   
            r_d = sqrt(((i_p - i)*dx).^2 + ((j_p - j)*dz).^2 + ((1-y)*dy).^2);
            
            if r_s > LFP_radius && r_d > LFP_radius
                continue
            end
            
            if r_s == 0 || r_d == 0
                continue
            end
            
            if r_s > LFP_radius
                I_p_s = 0; 
            else
                I_p_s = I_mat_PS(i_p,1,j_p,t_ind);
            end
            
            if r_d > LFP_radius
                I_p_d = 0; 
            else
                I_p_d = I_mat_PS(i_p,2,j_p,t_ind);
            end
            
            temp_sum_s = I_p_s/r_s + I_p_d/r_d;
            sum_eP = sum_eP + temp_sum_s;
        end
    end

    Ve_vec(t_ind) = gamma*sum_eP;
end

VoltVec = Ve_vec;

[~, LFP_vec] = convertSignal(VoltVec, sf);

Ve_vec = zeros(len_Ts, 1);
MUA_radius = 0.01;

for t_ind = 1:len_Ts
    
    sum_eP = 0;
    for i_p = 1:Nx
        for j_p = 1:Nz
            
%             if (abs(i_p - i)*dx > MUA_radius) || (abs(j_p - j)*dz > MUA_radius)
%                 continue
%             end
            
            r_s = sqrt(((i_p - i)*dx).^2 + ((j_p - j)*dz).^2 + ((0-y)*dy).^2);   
            r_d = sqrt(((i_p - i)*dx).^2 + ((j_p - j)*dz).^2 + ((1-y)*dy).^2);
            
            if r_s > MUA_radius && r_d > MUA_radius
                continue
            end
            
            if r_s == 0 || r_d == 0
                continue
            end
            
            if r_s > MUA_radius
                I_p_s = 0; 
            else
                I_p_s = I_mat(i_p,1,j_p,t_ind);
            end
            
            if r_d > MUA_radius
                I_p_d = 0; 
            else
                I_p_d = I_mat(i_p,2,j_p,t_ind);
            end
            
            temp_sum_s = I_p_s/r_s + I_p_d/r_d;
            sum_eP = sum_eP + temp_sum_s;
        end
    end

    Ve_vec(t_ind) = gamma*sum_eP;
end

MUA_VoltVec = Ve_vec;

[MUA_vec, ~] = convertSignal(MUA_VoltVec, sf);

end
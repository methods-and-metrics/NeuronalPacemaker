rng(1)

filetype = "BS_20by10_DGG_GABAton0-1_GABA_A0-1_test_31sec";
output_data = 0;
output_folder = "D:\BS_Model_Paper\Data\" + filetype + "\";
mkdir(char(output_folder));

% Constant Parameters
g_K_s = 0.077; % For dendrites: Persistent Potassium channel conductace in mS/cm^2 | Source: (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/
g_Na_p = 0.022; % For dendrites Persistent sodium channel conductace in mS/cm^2 | Source: (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/
g_K_d = 15; % For soma: potassium channel conductance in mS/cm^2 | See: https://www.mdpi.com/2073-4409/12/18/2229
g_Na = 30; % For soma: sodium channel conductance in mS/cm^2 | See: https://www.mdpi.com/2073-4409/12/18/2229
g_leak_d = 0.025; % For dendrite: leak conductance in mS/cm^2 | Source (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/
g_leak_s = 0.1; % For soma: sodium channel conductance in mS/cm^2 | See: https://www.mdpi.com/2073-4409/12/18/2229
g_Ca_max = 0.1;
g_K_Ca = 0; % will be deleted
g_C_Ca = 0; % will be deleted
p_s = 0.5; % Relative weight between soma and dendrite for intra-cellular conduction (A.u.)
g_c = 0.001; % Conductivity of soma-dendrite coupling in mS/cm^2
g_c_s = g_c/p_s; % scaled conductivity for soma
g_c_d = g_c/(1-p_s); % scaled conductivity for dendrite

D_GABA = 0.235;
g_GABA_ton = 0.5;
g_GABA_ton_scaling = 1;
minimal_g_GABA_A = 0.1; % Conductivity of GABA-A receptor in mS/cm^2 | Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2605954/#R43
g_GABA_B = 0.1; % Conductivity of GABA-B receptor in mS/cm^2 | Source (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC2605954/#R43

recovery_curr_scaling = 0.02;
% recovery_curr_scaling = 0.02;

lambda = 3; % In Hz, the spiking rate for glutamergic input
g_AMPA = 0.2; % Conductivity of AMPA channels in ms/cm^2

E_K = -85; % Reversal potential of potassium channel in mV | See (adapted): https://www.mdpi.com/2073-4409/12/18/2229
E_Na = 50; % Reversal potential of sodium channel in mV | See (adapted): https://www.mdpi.com/2073-4409/12/18/2229
E_leak = -70; % Reversal potential of leak current in mV | See (adapted): https://www.mdpi.com/2073-4409/12/18/2229
E_Ca = 75; % Delete
E_Cl = -70; % Reversal potential of chloride channels in mV | Source (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC2605954/#R43
E_GABA = -70; % Reversal potential for tonic inhibition, non outward rectifying
E_AMPA = 50; % Reversal potential of AMPA channel in mV

T = 32; % Temperature of the brain in Celsuis (A parameter for pacemaking) | Source: (adapted): https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/
C = 0.25; % Capicitance of dendrite in uF/cm^2 | Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/

V_offset = 70; % Voltage offset to calculate gating equations for pinsky model (this model shifts voltage so that resting potential is 0 mV)

neuron_density = 57000; % in units of neurons/mm^3 source: https://www.science.org/doi/10.1126/science.adk4858

% Changing Parameters
EF = 1;  % Binary on whether to include ephaptic effects  

% Number of grid points
Nx = 20;
Nz = 10;
Ny = 2;

% Dimensions of grid points in microns
dx = 50;
dz = 50;
dy = 250;

% num_synapses = floor(num_neurons*frac_synapse);

% Create create synaptic network
num_synapses = 5;
num_neurons = Nx*Nz;
synapse_inds = zeros(Nx, Nz, num_synapses);

GABA_alpha = 1;

for i = 1:Nx
    for j = 1:Nz
        for k = 1:num_synapses
            temp_ind = randi(num_neurons) - 1;
            temp_x_ind = floor(temp_ind/Nx);
            temp_z_ind = temp_ind - temp_x_ind*Nx;
            
            x_ind = temp_x_ind + 1;
            z_ind = temp_z_ind + 1;
            
            while (i == x_ind) && (j == z_ind)
                temp_ind = randi(num_neurons) - 1;
                temp_x_ind = floor(temp_ind/Nx);
                temp_z_ind = temp_ind - temp_x_ind*Nx;

                x_ind = temp_x_ind + 1;
                z_ind = temp_z_ind + 1;
            end
            
            synapse_inds(i, j, k) = temp_ind;
        end
    end
end

SF = neuron_density*((dx/1000)*(dz/1000)*(dy/1000)); % Adjust stacking factor to represent number of neurons within each grid point
p = 380;                                                               % In units of ohms*cm Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6312416/  
gamma = SF*p/(4*pi);                                                   % Scaling constant to convert current to extra-cellular potential | See https://www.mdpi.com/2073-4409/12/18/2229

soma_area = 400e-8; % This is assuming a soma diameter of around 20 microns, which would be considered a medium sized soma surface area for a pyramidal cell in the neocortex
dendrite_area = 36333e-8;  % In units of cm^2 Source: https://academic.oup.com/cercor/article/27/11/5398/4159219

T_dur = 10000; % Number of timesteps
dt = 0.1; % time step in milliseconds
sf = 1000*(1/dt); % sampling frequency

T_list = 0:dt:T_dur;

num_tpoints = length(T_list);

% Initialize matrices for calculations 
n_p_mat = zeros(Nx, Nz, num_tpoints);
n_mat = zeros(Nx, Nz, num_tpoints);
m_p_mat = zeros(Nx, Nz, num_tpoints);
h_mat = zeros(Nx, Nz, num_tpoints);
V_mat = ones(Nx, Ny, Nz, num_tpoints)*-70;
I_mat = zeros(Nx, Ny, Nz, num_tpoints);
I_mat_PS = zeros(Nx, Ny, Nz, num_tpoints);
Ve_mat = zeros(Nx, Ny, Nz, num_tpoints);
GABA_spike_time_mat = zeros(Nx, Nz);
spike_indicator_mat = zeros(Nx, Nz);
spike_timestamp_cells = cell(Nx, Nz);
Ca_mat = zeros(Nx, Nz, num_tpoints);
ex_Ca_vec = zeros(1,num_tpoints);
spike_train_mat = poissrnd(lambda*dt*1e-3, Nx, Nz, num_tpoints);
syn_spike_times = zeros(Nx, Nz);
GABA_tonic_vec = ones(1, num_tpoints)*0.1;
g_Ca_vec = ones(1,num_tpoints)*g_Ca_max;
GABA_tonic_mat = ones(Nx, Nz, num_tpoints);
recovery_curr_vec = zeros(1, num_tpoints);
f_m_temp_vec = zeros(1, num_tpoints);
fr_proxy_measure = zeros(Nx, Nz);

V_mat(:,:,:,1) = -65;

rand_T_var = 10*rand(Nx,Nz);
rand_inj_curr = 0.1*rand(Nx, Nz);

% Start each cell at different times to ensure subthreshold oscillations
% are out of phase for each cell
rand_start_times = floor(rand(Nx, Nz)*10000);

% Convert from microns to cm, for electric potential calculation
cm_dx = dx*1e-4;
cm_dy = dy*1e-4;
cm_dz = dz*1e-4;

I_inj = 0.8; % Baseline injected current into dendrites to create intrinsic subthreshold oscillations
soma_current = 0.05; % Injected current into soma to simulate excitation (and the propensity for intrinsic firing)

tonic_time_window = 2; % In seconds

Ca_spike_avail = 1;
GABA_spike_dist = 0; % Initialize spike time distance variable
Glu_spike_dist = 0;
last_spike_time_Ca = 0;
tonic_spike_lag = 1000; % In milliseconds

temp_oscillations = 0.5 + (1+sin((T_list/1000)*2*pi*(1/2)))/2;

basal_firing = 0;

GABA_AP_factor = 5.7e-4;

for t = 1:(num_tpoints-1)
    
    disp(t)
    temp_Ca_conc = 0;
    temp_GABA_conc = GABA_tonic_vec(t);
    f_m_temp = 0;
    for i = 1:Nx
        for j = 1:Nz
            
            if rand_start_times(i,j) > t
                continue
            end
            
            if t <= 10000
                g_AMPA = 0.1;
                g_GABA_ton_temp = GABA_alpha*temp_GABA_conc;
                GABA_AP_factor = 5.7e-4;
                firing_rate_time_factor = (1/1000);
                soma_current = 0.06;
                basal_firing = 0.1;
                minimal_g_GABA_A = 0.3;
            elseif t <= 310000
                g_AMPA = 0;
%                 g_GABA_ton_temp = (GABA_alpha+(0.03*(t-10000)/10000))*temp_GABA_conc;
                g_GABA_ton_temp = GABA_alpha*temp_GABA_conc;
                GABA_AP_factor = 5.7e-4;
                firing_rate_time_factor = (1/1000);
                soma_current = 0.06;
                basal_firing = 0.1;
                minimal_g_GABA_A = 0.3;
            else
                g_AMPA = 0;
%                 g_GABA_ton_temp = (GABA_alpha+(0.03*(t-10000)/10000))*temp_GABA_conc;
                g_GABA_ton_temp = GABA_alpha*(10)*temp_GABA_conc;
                GABA_AP_factor = 5.7e-4;
                firing_rate_time_factor = (1/1000);
                soma_current = 0.06;
                basal_firing = 0;
                minimal_g_GABA_A = 0.1;
            end
%             else
%                 g_AMPA = 0.1;
%                 g_GABA_ton_temp = GABA_alpha*(0.1)*temp_GABA_conc;
%                 soma_current = 0.06;
%                 basil_firing = 0.1;
%                 minimal_g_GABA_A = 0;
%             end
            
            n_p = n_p_mat(i,j,t);
            n = n_mat(i,j,t);
            m_p = m_p_mat(i,j,t);
            h = h_mat(i,j,t);
            V_s = V_mat(i,1,j,t);
            V_d = V_mat(i,2,j,t);
            Ve_s = Ve_mat(i,1,j,t);
            Ve_d = Ve_mat(i,2,j,t);
            Ca_conc = Ca_mat(i,j,t);
            ex_Ca_conc = ex_Ca_vec(t);
            spike_indicator = spike_train_mat(i,j,t);
            temp_g_Ca = g_Ca_vec(t);
            rc_value = recovery_curr_vec(t);
            temp_firing_rate = fr_proxy_measure(i,j);
            
            if spike_indicator > 0
                syn_spike_times(i,j) = t;
            end
            
            last_spike_time = GABA_spike_time_mat(i, j);
            spike_indicator = spike_indicator_mat(i, j);
            temp_spike_timestamps = spike_timestamp_cells{i, j};
            
            if last_spike_time > 0
                GABA_spike_dist = (t-last_spike_time)*dt;
            else
                GABA_spike_dist = 0;
            end
%             disp(GABA_spike_dist)
            
            Glu_spike_time = syn_spike_times(i,j);
            
            if Glu_spike_time > 0
                Glu_spike_dist = (t-Glu_spike_time)*dt;
            else
                Glu_spike_dist = 0;
            end
            
            if t > 1
                dVe_s = Ve_s - Ve_mat(i,1,j,t-1);
                dVe_d = Ve_d - Ve_mat(i,2,j,t-1);
            else
                dVe_s = 0;
                dVe_d = 0;
            end

            if ~isempty(temp_spike_timestamps)
                oldest_spiketime = temp_spike_timestamps(1);
                if (t - oldest_spiketime) > (tonic_time_window*(1000/dt))
                    if length(temp_spike_timestamps) == 1
                        temp_spike_timestamps = [];
                    else
                        temp_spike_timestamps = temp_spike_timestamps(2:end);
                        spike_timestamp_cells{i, j} = temp_spike_timestamps;
                    end
                end
            end
            
            g_GABA_A = minimal_g_GABA_A;
            
            Ca_weight = max(ex_Ca_conc,0);
            
            if ex_Ca_conc < 0.1
                Ca_spike_avail = 0;
            else
                g_Ca = temp_g_Ca + dt*(g_Ca_max - temp_g_Ca);
            end
            
            if Ca_spike_avail == 0
                g_Ca = temp_g_Ca + dt*(-temp_g_Ca)/100;
                if ex_Ca_conc > 0.9
                    Ca_spike_avail = 1;
                end
            end
            
            g_Ca_vec(t+1) = g_Ca;
            
            synapse_spike_vec = zeros(num_synapses, 1);
            num_spikes = length(temp_spike_timestamps);
%             if num_spikes > 0
%                 for sp = 1:num_spikes
% %                     conv_spike_dist = dt*(t - temp_spike_timestamps(sp));
% %                     conv_spike_dist_lag = max(0, conv_spike_dist-tonic_spike_lag);
% %                     local_ton_GABA_conc = local_ton_GABA_conc + g_GABA_ton_scaling*(1 - exp(-conv_spike_dist_lag/200))*exp(-conv_spike_dist_lag/1000);
%                 end
%             end
            new_firing_rate = temp_firing_rate + dt*(firing_rate_time_factor*(num_spikes/tonic_time_window - temp_firing_rate));

            f_m_temp = f_m_temp +  new_firing_rate;

            if num_spikes > 0
                recovery_curr = rc_value + dt*(num_spikes*recovery_curr_scaling - rc_value)/100;
            else
                recovery_curr = rc_value + dt*(-rc_value)/100;
            end
            
            
%             g_GABA_ton_temp = temp_oscillations(t);
            
            for syn_num = 1:num_synapses
                temp_ind = synapse_inds(i, j, syn_num);
                
                temp_z_ind = floor(temp_ind/Nx);
                temp_x_ind = temp_ind - temp_z_ind*Nx;

                x_ind = temp_x_ind + 1;
                z_ind = temp_z_ind + 1;
                
                synapse_spike_vec(syn_num) = dt*(t - GABA_spike_time_mat(x_ind, z_ind));
            end

            g_Ca = g_Ca_max;
            [V_s_new, n_new, h_new, I_s_trans, g_GABA_A_bar, I_s_trans_PS] = soma_neuron(-0.5 + soma_current*rand(1), V_s, V_d, Ve_s, dVe_s, n, h, dt, ...
                g_K_d, g_Na, g_leak_s, g_c_s, g_GABA_A, g_GABA_ton_temp, g_AMPA, E_K, E_Na, E_Cl, E_GABA, E_AMPA, E_leak, 3, V_offset, GABA_spike_dist, Glu_spike_dist, synapse_spike_vec);

            [V_d_new, n_p_new, m_p_new, I_d_trans, I_syn_d, Ca_conc_new, I_d_trans_PS] = dendrite_neuron(I_inj, V_d, V_d, Ve_d, dVe_d, n_p, m_p, dt, g_Ca, g_K_Ca, ...
                g_C_Ca, g_leak_d, g_K_s, g_Na_p, g_c_d, g_GABA_B, g_AMPA, E_Ca, E_K, E_Na, E_AMPA, E_leak, C, T - rand_T_var(i,j), GABA_spike_dist, Glu_spike_dist, Ca_conc, [], Ca_spike_avail);
            
            if i == 20 && j == 5
                disp(V_s_new)
                disp(V_d_new)
            end
            
            if spike_indicator == 1
                if V_s_new < 0
                    spike_indicator_mat(i,j) = 0;
                    spike_timestamp_cells{i, j} = [temp_spike_timestamps, t];
                end
            end
            
            if V_s_new > 0
                GABA_spike_time_mat(i,j) = t;
                spike_indicator_mat(i,j) = 1;
%                 temp_GABA_conc = temp_GABA_conc + 1/num_neurons;
            end
            
%             if i == 1 && j == 1
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(Nx, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, Nz, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             elseif i == 1 && j == Nz
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(Nx, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, 1, t))/dz^2);
%             elseif i == Nx && j == Nz
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, 1, t))/dz^2);
%             elseif i == Nx && j == 1
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i-1, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, Nz, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             elseif i == 1 && j < Nz
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(Nx, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             elseif i == Nx && j < Nz
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i-1, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             elseif i < Nx && j == 1
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i-1, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, Nz, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             elseif i < Nx && j == Nz
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i-1, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, 1, t))/dz^2);
%             else
%                 ex_GABA_diff = D_GABA*((GABA_tonic_mat(i-1, j, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i+1, j, t))/dx^2 + ...
%                     (GABA_tonic_mat(i, j-1, t) - 2*GABA_tonic_mat(i, j, t) + GABA_tonic_mat(i, j+1, t))/dz^2);
%             end

            n_p_mat(i,j,t+1) = n_p_new;
            n_mat(i,j,t+1) = n_new;
            m_p_mat(i,j,t+1) = m_p_new;
            h_mat(i,j,t+1) = h_new;
            V_mat(i,1,j,t+1) = V_s_new;
            V_mat(i,2,j,t+1) = V_d_new;
            I_mat(i,1,j,t+1) = I_s_trans*1e-3*soma_area;
            I_mat(i,2,j,t+1) = I_d_trans*1e-3*dendrite_area;
            I_mat_PS(i,1,j,t+1) = I_s_trans_PS*1e-3*soma_area;
            I_mat_PS(i,2,j,t+1) = I_d_trans_PS*1e-3*dendrite_area;
            Ca_mat(i,j,t+1) = Ca_conc_new;
            recovery_curr_vec(t+1) = recovery_curr;
            fr_proxy_measure(i,j) = new_firing_rate;
            
            temp_Ca_conc = temp_Ca_conc + (1-Ca_conc_new);
            
            sum_eP_s = 0;
            sum_eP_d = 0;
            if EF == 1
                if t >= 10000
                    for i_p = 1:Nx
                        for j_p = 1:Nz
                            if j_p == j && i_p == i
                                continue
                            end

                            r_same = sqrt(((i_p - i)*cm_dx).^2 + ((j_p - j)*cm_dz).^2);
                            r_opp = sqrt(((i_p - i)*cm_dx).^2 + ((j_p - j)*cm_dz).^2 + cm_dy.^2);

                            I_s_p = I_mat(i_p, 1, j_p, t);
                            I_d_p = I_mat(i_p, 2, j_p, t);

                            temp_sum_s = I_s_p/r_same + I_d_p/r_opp;
                            temp_sum_d = I_s_p/r_opp + I_d_p/r_same;
                            sum_eP_s = sum_eP_s + temp_sum_s;
                            sum_eP_d = sum_eP_d + temp_sum_d;
                        end
                    end
                end
            end
            
            Ve_mat(i,1,j,t+1) = gamma*sum_eP_s;
            Ve_mat(i,2,j,t+1) = gamma*sum_eP_d;
        end
    end
    
    ex_Ca_vec(t+1) = temp_Ca_conc/num_neurons;
    GABA_tonic_vec(t+1) = temp_GABA_conc + dt*(GABA_AP_factor*(f_m_temp/num_neurons + basal_firing) - 0.004*(GABA_tonic_vec(t)-0.1));
    f_m_temp_vec(t+1) = f_m_temp;
end

[MUA_vec, LFP_vec, VoltVec] = getMUA_soma_neuron(I_mat, I_mat_PS, num_tpoints, Nx, Ny, Nz, cm_dx, cm_dy, cm_dz, Nx/2, Nz/2, 1, gamma, sf);

figure(1)
clf
temp_array1 = zeros(Nx, num_tpoints);
temp_array1(:) = V_mat(:,1,floor(Nz/2),:);
mesh(temp_array1)
view(2)
if output_data
    temp_fn = "Vm_mat_soma.fig";
    savefig(char(output_folder+temp_fn));
end

figure(2)
clf
temp_array1 = zeros(Nx, num_tpoints);
temp_array1(:) = V_mat(:,2,floor(Nz/2),:);
mesh(temp_array1)
view(2);
if output_data
    temp_fn = "Vm_mat_dendrite.fig";
    savefig(char(output_folder+temp_fn));
end
    
figure(3)
clf
temp_array1 = zeros(Nx, num_tpoints);
temp_array1(:) = I_mat(:,1,floor(Nz/2),:);
mesh(temp_array1)
view(2);

figure(4)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = Ve_mat(floor(Nx/2),1,floor(Nz/2),:);
plot(T_list, temp_array1)

figure(5)
clf
plot(T_list, MUA_vec, "black")
% ylim([-0.05, 0.05])
if output_data
    temp_fn = "MUA_vec.fig";
    savefig(char(output_folder+temp_fn));
end

figure(6)
clf
plot(T_list, LFP_vec, "black")
ylim([-3, 3])
if output_data
    temp_fn = "LFP_vec.fig";
    savefig(char(output_folder+temp_fn));
end

figure(7)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = V_mat(floor(Nx/2),1,floor(Nz/2),:);
plot(T_list, temp_array1)

figure(9)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = I_mat(floor(Nx/2),2,floor(Nz/2),:);
plot(T_list, temp_array1)
if output_data
    temp_fn = "I_vec_center_neuron.fig";
    savefig(char(output_folder+temp_fn));
end

figure(10)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = V_mat(floor(Nx/2),2,floor(Nz/2),:);
plot(T_list, temp_array1)

figure(11)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = Ca_mat(floor(Nx/2),floor(Nz/2),:);
plot(T_list, temp_array1)

figure(12)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = ex_Ca_vec;
plot(T_list, temp_array1)
if output_data
    temp_fn = "ex_Ca_vec.fig";
    savefig(char(output_folder+temp_fn));
end

figure(13)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = GABA_tonic_vec;
plot(T_list, temp_array1)
if output_data
    temp_fn = "GABA_tonic_vec.fig";
    savefig(char(output_folder+temp_fn));
end

figure(14)
clf
temp_array1 = zeros(1, num_tpoints);
temp_array1(:) = f_m_temp_vec;
plot(T_list, temp_array1/num_neurons)
if output_data
    temp_fn = "avg_firing_vec.fig";
    savefig(char(output_folder+temp_fn));
end
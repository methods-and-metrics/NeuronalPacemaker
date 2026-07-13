function [V_new, n_new, h_new, I_trans, g_GABA_A_bar, I_trans_PS] = soma_neuron(I_inj, V_s, V_d, Ve, dVe, n, h, dt, ...
    g_K_d, g_Na, g_leak, g_c, g_GABA_A, g_GABA_ton, g_AMPA, E_K, E_Na, E_Cl, E_GABA, E_AMPA, E_leak, C, V_offset, GABA_spike_dist, Glu_spike_dist, synapse_spike_vec)

    g_Na_bar = g_Na*(m_ss(V_s + V_offset)^2)*h;
    g_K_d_bar = g_K_d*n;
    g_GABA_A_bar = g_GABA_A*(1 - exp(-GABA_spike_dist/2))*exp(-GABA_spike_dist/50);
    g_GABA_ton_bar = g_GABA_ton;
    
    num_synapses = length(synapse_spike_vec);
    g_AMPA_bar = 0;
    for i = 1:num_synapses
        g_AMPA_bar = g_AMPA_bar + (g_AMPA*(1 - exp(-synapse_spike_vec(i)/0.3))*exp(-synapse_spike_vec(i)/3))/num_synapses;
    end
    g_AMPA_bar = g_AMPA_bar + g_AMPA*(1 - exp(-Glu_spike_dist/0.3))*exp(-Glu_spike_dist/3);

    I_K_d = g_K_d_bar*(V_s - E_K);
    I_Na = g_Na_bar*(V_s - E_Na);
    I_leak = g_leak*(V_s - E_leak);
    I_ds = g_c*(V_d - V_s);
    I_GABA_A = g_GABA_A_bar*(V_s - E_Cl);
    I_GABA_ton = g_GABA_ton_bar*(V_s - E_GABA);
    I_AMPA = g_AMPA_bar*(V_s - E_AMPA);
    
    I_ion = I_K_d + I_Na + I_leak;
    I_syn = I_GABA_A + I_GABA_ton + I_AMPA;
    I_eph_res = (g_K_d_bar + g_Na_bar + g_leak)*Ve;
    
    V_new = V_s - dVe + (dt*(-I_ion - I_eph_res - I_syn + I_inj + I_ds))/C;

    h_new = h + dt*(h_ss(V_s+V_offset) - h)/h_tau(V_s+V_offset);
    n_new = n + dt*(n_ss(V_s+V_offset) - n)/n_tau(V_s+V_offset);

    I_trans = I_ion - I_inj + I_eph_res + I_syn;
    I_trans_PS = I_syn + I_ion;
end

% Try new Sodium and potassium channel conductances... in accordance to https://pmc.ncbi.nlm.nih.gov/articles/PMC2605954/
% Only keep the persistent potassium and sodium channels from https://pmc.ncbi.nlm.nih.gov/articles/PMC1157807/pdf/jphysiol00327-0073.pdf

function hss = h_ss(Vs)
    hss = alphaH(Vs)/(alphaH(Vs) + betaH(Vs));
end

function htau = h_tau(Vs)
    htau = 1/(alphaH(Vs) + betaH(Vs));
end

function nss = n_ss(Vs)
    nss = alphaN(Vs)/(alphaN(Vs) + betaN(Vs));
end

function ntau = n_tau(Vs)
    ntau = 1/(alphaN(Vs) + betaN(Vs));
end

function mss = m_ss(Vs)
    mss = alphaM(Vs)/(alphaM(Vs) + betaM(Vs));
end

% % Parameters for steady state functions

function aH = alphaH(Vs)
    aH = 0.128*exp((17-Vs)/18);
end

function bH = betaH(Vs)
    bH = 4/(1 + exp((40-Vs)/5));
end

function aN = alphaN(Vs)
    aN = (0.016*(35.1-Vs))/(exp((35.1 - Vs)/5) - 1);
end

function bN = betaN(Vs)
    bN = 0.25*exp(0.5 - 0.025*Vs);
end

function aM = alphaM(Vs)
    aM = (0.32*(13.1-Vs))/(exp((13.1-Vs)/4) - 1);
end

function bM = betaM(Vs)
    bM = (0.28*(Vs-40.1))/(exp((Vs-40.1)/5) - 1);
end

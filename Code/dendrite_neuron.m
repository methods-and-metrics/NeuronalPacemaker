function [V_new, n_p_new, m_p_new, I_trans, I_Ca, Ca_conc, I_trans_PS] = dendrite_neuron(I_inj, V_d, V_s, Ve, dVe, n_p, m_p, dt, g_Ca, g_K_Ca, ...
    g_C_Ca, g_leak, g_K_s, g_Na_p, g_c, g_GABA_B, g_AMPA, E_Ca, E_K, E_Na, E_AMPA, E_leak, C, T, GABA_spike_dist, Glu_spike_dist, Ca_conc, synapse_spike_vec, Ca_spike_avail)
    
    Ca_alpha_function = exp(1)*(GABA_spike_dist/10)*exp(-GABA_spike_dist/10);
    
    g_AMPA_bar = 0;
    
    num_synapses = length(synapse_spike_vec);
    
    for i = 1:num_synapses
        g_AMPA_bar = g_AMPA_bar + (g_AMPA*(1 - exp(-synapse_spike_vec(i)/0.3))*exp(-synapse_spike_vec(i)/3))/num_synapses;
    end

    g_K_s_bar = g_K_s*n_p;
    g_Na_p_bar = g_Na_p*m_p;
    if GABA_spike_dist < 200
        g_GABA_B_bar = g_GABA_B*(1 - exp(-GABA_spike_dist/50));
    else
        g_GABA_B_bar = exp(-GABA_spike_dist/50);
    end

    I_Ca = g_Ca*Ca_alpha_function*(V_d - E_Ca);
    I_leak = g_leak*(V_d - E_leak);
    I_ds = g_c*(V_d - V_s);

    I_K_s = g_K_s_bar*(V_d - E_K);
    I_Na_p = g_Na_p_bar*(V_d - E_Na);
    
    I_GABA_B = g_GABA_B_bar*(V_d - E_K);
    I_AMPA = g_AMPA_bar*(V_d - E_AMPA);

    I_ion = I_leak + I_K_s + I_Na_p + I_Ca;
    I_syn = I_GABA_B + I_AMPA;
    I_eph_res = (g_K_s_bar + g_Na_p_bar + g_leak)*Ve;
    
    V_new = V_d - dVe + dt*(-I_ion + I_inj - I_eph_res - I_syn - I_ds)/C;
    
    n_p_inf = n_p_inf_func(V_d);
    tau_K_s = tau_K_s_func(V_d, T);
    m_p_inf = m_p_inf_func(V_d);
    tau_Na_p = tau_Na_p_func();

    n_p_new = n_p + dt*((n_p_inf - n_p)/tau_K_s);
    m_p_new = m_p + dt*((m_p_inf - m_p)/tau_Na_p);
    
    Ca_conc = Ca_conc + dt*(-0.001*I_Ca - (1/1000)*Ca_conc);
    
    I_trans = I_ion - I_inj + I_syn + I_eph_res;
    I_trans_PS = I_syn + I_ion;
    
end

function ninf = n_p_inf_func(V)
    ninf = 1/(1 + exp(-(V+35)/10));
end

function tauKs = tau_K_s_func(V, T)
    tauKs = (1000/(3.3*(exp((V+35)/40) + exp(-(V+35)/20))))*(1/(3^((T-22)/10)));
end

function mpinf = m_p_inf_func(V)
    mpinf = 1/(1 + exp(-(V+40)/5));
end

function tauNaP = tau_Na_p_func()
    tauNaP = 5;
end
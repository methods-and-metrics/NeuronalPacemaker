function [V_new, n_p_new, m_p_new, I_trans, I_Ca, Ca_conc, I_trans_PS] = dendrite_neuron(I_inj, V_d, V_s, Ve, dVe, n_p, m_p, dt, g_Ca, g_K_Ca, ...
    g_C_Ca, g_leak, g_K_s, g_Na_p, g_c, g_GABA_B, g_AMPA, E_Ca, E_K, E_Na, E_AMPA, E_leak, C, T, GABA_spike_dist, Glu_spike_dist, Ca_conc, synapse_spike_vec, Ca_spike_avail)

    S_depr = 1;
    
%     Ca_spike_dur = 50; % in milliseconds
    
%     if (GABA_spike_dist > 0 && GABA_spike_dist < Ca_spike_dur)
% %         Ca_indicator = 2/(exp((GABA_spike_dist-25)/8) + exp(-(GABA_spike_dist-25)/8));
% %         Ca_conc = Ca_conc + dt*(Ca_spike_level/Ca_spike_dur);
%         
%         Ca_alpha_function
%     else
%         Ca_indicator = 0;
%     end
    
    Ca_alpha_function = exp(1)*(GABA_spike_dist/10)*exp(-GABA_spike_dist/10);
    
    g_AMPA_bar = 0;
    
    num_synapses = length(synapse_spike_vec);
    
    for i = 1:num_synapses
%         g_AMPA_bar = g_AMPA_bar + (1/10)*(g_AMPA*(1 - exp(-synapse_spike_vec(i)/0.3))*exp(-synapse_spike_vec(i)/3))/num_synapses;
        g_AMPA_bar = g_AMPA_bar + (g_AMPA*(1 - exp(-synapse_spike_vec(i)/0.3))*exp(-synapse_spike_vec(i)/3))/num_synapses;
    end

%     g_AMPA_bar = g_AMPA_bar + g_AMPA*(1 - exp(-Glu_spike_dist/0.3))*exp(-Glu_spike_dist/3);

    g_K_s_bar = g_K_s*n_p;
    g_Na_p_bar = g_Na_p*m_p;
    if GABA_spike_dist < 200
        g_GABA_B_bar = g_GABA_B*(1 - exp(-GABA_spike_dist/50));%*S_depr*(1-exp(-GABA_spike_dist/2))*exp(-GABA_spike_dist/100);
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
%     I_NMDA = g_NMDA_bar*(V_d - E_NMDA);

    I_ion = I_leak + I_K_s + I_Na_p + I_Ca;
    I_syn = I_GABA_B + I_AMPA;
    I_eph_res = (g_K_s_bar + g_Na_p_bar + g_leak)*Ve;
    
    V_new = V_d - dVe + dt*(-I_ion + I_inj - I_eph_res - I_syn - I_ds)/C;
    
%     alpha_s = alpha_s_func(V_d);
%     beta_s = beta_s_func(V_d);
%     alpha_r = alpha_r_func(V_d);
%     beta_r = beta_r_func(V_d, alpha_r);
%     alpha_q = alpha_q_func();
%     beta_q = beta_q_func();
%     alpha_c = alpha_c_func(V_d);
%     beta_c = beta_c_func(V_d, alpha_c);

    n_p_inf = n_p_inf_func(V_d);
    tau_K_s = tau_K_s_func(V_d, T);
    m_p_inf = m_p_inf_func(V_d);
    tau_Na_p = tau_Na_p_func();
    
%     s_new = s + dt*(alpha_s*(1-s) - beta_s*s);
%     r_new = r + dt*(alpha_r*(1-r) - beta_r*r);
%     q_new = q + dt*(alpha_q*(1-q) - beta_q*q);
%     c_new = c + dt*(alpha_c*(1-c) - beta_c*c);

    n_p_new = n_p + dt*((n_p_inf - n_p)/tau_K_s);
    m_p_new = m_p + dt*((m_p_inf - m_p)/tau_Na_p);
    
    Ca_conc = Ca_conc + dt*(-0.001*I_Ca - (1/1000)*Ca_conc);
    
    I_trans = I_ion - I_inj + I_syn + I_eph_res;
    I_trans_PS = I_syn + I_ion;
    
end

function alphas = alpha_s_func(V)
    alphas = 1.6/(1 + exp(-0.072*(V-5)));
end

function betas = beta_s_func(V)
    betas = (0.02*(V + 13.9))/(exp((V + 8.9)/5)-1);
end

function alphar = alpha_r_func(V)
    if V <= -65
        alphar = 0.005;
    else
        alphar = exp((-V + 65)/20)/200;
    end
end

function betar = beta_r_func(V, ar)
    if V <= -65
        betar = 0;
    else
        betar = 0.005 - ar;
    end
end

function alphaq = alpha_q_func()
    alphaq = 0.01;
end

function betaq = beta_q_func()
    betaq = 0.001;
end

function alphac = alpha_c_func(V)
    if V <= -15
        alphac = exp((V + 55)/11 - (V + 58.5)/27)/19;
    else
        alphac = 2*exp(-(V + 58.5)/27);
    end
end

function betac = beta_c_func(V, ac)
    if V <= -15
        betac = 2*exp(-(V + 58.5)/27) - ac;
    else
        betac = 0;
    end
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
"""
DA_ODE
--------------------
Fonction biophysique de base du modèle de Yu SANS contrôleur 
homéostatique. Intègre les canaux de type N (Cav2.2) et T (Cav3.1) 
en plus des canaux du modèle original. Les conductances sont des 
PARAMÈTRES FIXES — elles ne varient pas au cours du temps.

Utilisée principalement dans le notebook pour :
  - La calibration de gCaT (contribution à 10% du flux calcique)
  - Le calcul de la cible calcique (Ca_tgt) de chaque neurone
  - La vérification du comportement de base avant intégration 
    du contrôleur
"""


function DA_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1](t)  
    gNa   = p[2]     
    gCaL  = p[3]     
    gKd   = p[4]     
    gKA   = p[5]     
    gKERG = p[6]     
    gKSK  = p[7]     
    gH    = p[8]     
    gLNS  = p[9]     
    gLCa  = p[10]    
    gNtype= p[11]
    gTtype= p[12]

    # State variables
    V    = u[1]      
    m    = u[2]      
    h    = u[3]     
    hs   = u[4]      
    l    = u[5]      
    n    = u[6]      
    p_ka = u[7]      
    q1   = u[8]      
    q2   = u[9]     
    o    = u[10]     
    i    = u[11]     
    mH   = u[12]     
    Ca   = u[13]     
    mN, hN, mT, hT = u[14:17]
    
    # Calcium-dependent SK current and pump
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end

    # Membrane potential dynamics
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gCaL * l * (V - VCa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        + 100 * Iapp / (pi * d * L)
    )

    # Gating variable dynamics
    du[2]  = (1 / tau_m(V))  * (m_inf(V)  - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)  - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V) - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf(V)  - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)  - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)  - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V) - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V) - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[14] = (mN_inf(V) - mN) / tau_mN(V)
    du[15] = (hN_inf(V) - hN) / tau_hN(V)    
    du[16] = (mT_inf(V) - mT) / tau_mT(V)
    du[17] = (hT_inf(V) - hT) / tau_hT(V)
    
    # Calcium dynamics
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa)) / (F * d * 0.1)
end

"""
DA_homeo_2014_ODE_Ntype_Ttype_Ltypeblock_LeakC
------------------------------------------
Fonction principale pour simuler le BLOCAGE DES CANAUX L-TYPE (Cav1.3)
avec contrôleur homéostatique intégré.

Caractéristiques :
  - Canaux N-type et T-type actifs comme sources calciques alternatives
  - Conductance de fuite calcique (gLCa) FIXE — découplée du 
    contrôleur homéostatique 
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  
"""

function DA_homeo_2014_ODE_Ntype_Ttype_Ltypeblock_LeakC(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]    
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa, tNCa, tTCA = p[4:14]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances 
    g_max_main = p[15][1]
    g_max_cal = p[15][2]  
   
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = clamp.(u[14:22], 0.0, g_max_main)
    gNtype, gTtype = clamp.(u[32:33], 0.0, g_max_cal)
    
    ms_main = max.(0.0, u[23:31])
    ms_cal  = max.(0.0, u[34:35])
    
    mN, hN, mT, hT = u[36:39]

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    

    gCaL_regulee = u[15] 
    gCaL_eff = (t < 30000.0) ? gCaL_regulee : 1e-18
    tCaL_eff = (t < 30000.0) ? p[5] : 1e16

    # --- 4. DYNAMIQUE DU POTENTIEL MEMBRANAIRE ---
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gCaL_eff * l * (V - VCa) 
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        + 100 * Iapp / (pi * d * L)
    )

    # --- 5. DYNAMIQUE DES PORTES  ---
    du[2]  = (1 / tau_m(V))  * (m_inf(V)  - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)  - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V) - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf(V)  - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)  - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)  - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V) - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V) - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[36] = (mN_inf(V) - mN) / tau_mN(V)
    du[37] = (hN_inf(V) - hN) / tau_hN(V)    
    du[38] = (mT_inf(V) - mT) / tau_mT(V)
    du[39] = (hT_inf(V) - hT) / tau_hT(V)

    # --- 6. DYNAMIQUE DU CALCIUM  ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa)+ gTtype * mT * hT * (V - VCa))/ (F * d * 0.1)

    for k in 0:7 
        du[14+k] = (1/tau_g) * (u[23+k] - u[14+k]) 
    end
    
    du[22] = 0.0                       
    du[32] = (1/tau_g) * (u[34] - u[32])
    du[33] = (1/tau_g) * (u[35] - u[33])
    
    error = (Ca_tgt - Ca) 
    ts_main = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tLNS]
    for k in 0:7
        if (u[14+k] >= g_max_main[k+1] && error > 0) || (u[23+k] <= 0.0 && error < 0)
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts_main[k+1]) * error
        end
    end
    
    du[31] = 0.0
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if (u[32+k] >= g_max_cal[k+1] && error > 0) || (u[34+k] <= 0.0 && error < 0)
            du[34+k] = 0.0
        else
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end
  
end

"""
DA_homeo_2014_AType
---------------------
Fonction pour simuler le BLOCAGE DES CANAUX A-TYPE (Kv4)
avec contrôleur homéostatique intégré.

Caractéristiques :
  - Identique à DA_homeo_2014_ODE_Ntype_Ttype_Ltypeblock_LeakC
    dans sa structure générale
  - Blocage instantané de gKA à t = 30000 ms 
  - Tous les autres canaux, y compris gCaL, restent actifs 
    et régulés normalement
  - gLCa fixe, gLNS régulé normalement
"""


function DA_homeo_2014_AType(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]    
    tau_g      = p[3]   
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa, tNCa, tTCA = p[4:14]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances 
    g_max_main = p[15][1]
    g_max_cal = p[15][2]  
   
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = clamp.(u[14:22], 0.0, g_max_main)
    gNtype, gTtype = clamp.(u[32:33], 0.0, g_max_cal)
    
    ms_main = max.(0.0, u[23:31])
    ms_cal  = max.(0.0, u[34:35])
    
    mN, hN, mT, hT = u[36:39]

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
  
    gKA_regulee = u[17] 
    gKA_eff = (t < 30000.0) ? gKA_regulee : 1e-18
    tKA_eff = (t < 30000.0) ? p[7] : 1e16

    # --- 4. DYNAMIQUE DU POTENTIEL MEMBRANAIRE ---
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gCaL * l * (V - VCa) 
        - gKd * n^3 * (V - VK)
        - gKA_eff * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        + 100 * Iapp / (pi * d * L)
    )

    # --- 5. DYNAMIQUE DES PORTES  ---
    du[2]  = (1 / tau_m(V))  * (m_inf(V)  - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)  - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V) - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf(V)  - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)  - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)  - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V) - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V) - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[36] = (mN_inf(V) - mN) / tau_mN(V)
    du[37] = (hN_inf(V) - hN) / tau_hN(V)    
    du[38] = (mT_inf(V) - mT) / tau_mT(V)
    du[39] = (hT_inf(V) - hT) / tau_hT(V)

    # --- 6. DYNAMIQUE DU CALCIUM ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL* l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa)+ gTtype * mT * hT * (V - VCa))/ (F * d * 0.1)

    for k in 0:7 
        du[14+k] = (1/tau_g) * (u[23+k] - u[14+k]) 
    end
    
    du[22] = 0.0                           

    du[32] = (1/tau_g) * (u[34] - u[32])
    du[33] = (1/tau_g) * (u[35] - u[33])
    
    error = (Ca_tgt - Ca) 
    ts_main = [t_Na, tCaL, tKd, tKA_eff, tKERG, tKSK, tH, tLNS]
    for k in 0:7
        if (u[14+k] >= g_max_main[k+1] && error > 0) || (u[23+k] <= 0.0 && error < 0)
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts_main[k+1]) * error
        end
    end
    
    du[31] = 0.0
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if (u[32+k] >= g_max_cal[k+1] && error > 0) || (u[34+k] <= 0.0 && error < 0)
            du[34+k] = 0.0
        else
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end
  
end

"""
DA_homeo_2014_HType
---------------------
Fonction pour simuler le BLOCAGE DES CANAUX HCN (Ih)
avec contrôleur homéostatique intégré.

Caractéristiques :
  - Identique à DA_homeo_2014_ODE_Ntype_Ttype_Ltypeblock_LeakC
    dans sa structure générale
  - Blocage instantané de gH à t = 30000 ms 
  - Tous les autres canaux, y compris gCaL, restent actifs 
    et régulés normalement
  - gLCa fixe, gLNS régulé normalement
"""
function DA_homeo_2014_HType(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]   
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa, tNCa, tTCA = p[4:14]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13])
    
    # Conductances 
    g_max_main = p[15][1]
    g_max_cal = p[15][2]  
   
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = clamp.(u[14:22], 0.0, g_max_main)
    gNtype, gTtype = clamp.(u[32:33], 0.0, g_max_cal)
    
    ms_main = max.(0.0, u[23:31])
    ms_cal  = max.(0.0, u[34:35])
    
    mN, hN, mT, hT = u[36:39]

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    
    gH_regulee = u[20] 
    gH_eff = (t < 30000.0) ? gH_regulee : 1e-18
    tH_eff = (t < 30000.0) ? p[10] : 1e16

    # --- 4. DYNAMIQUE DU POTENTIEL MEMBRANAIRE ---
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gCaL * l * (V - VCa) 
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH_eff * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        + 100 * Iapp / (pi * d * L)
    )

    # --- 5. DYNAMIQUE DES PORTES ---
    du[2]  = (1 / tau_m(V))  * (m_inf(V)  - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)  - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V) - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf(V)  - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)  - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)  - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V) - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V) - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[36] = (mN_inf(V) - mN) / tau_mN(V)
    du[37] = (hN_inf(V) - hN) / tau_hN(V)    
    du[38] = (mT_inf(V) - mT) / tau_mT(V)
    du[39] = (hT_inf(V) - hT) / tau_hT(V)

    # --- 6. DYNAMIQUE DU CALCIUM ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL* l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa)+ gTtype * mT * hT * (V - VCa))/ (F * d * 0.1)

    for k in 0:7 
        du[14+k] = (1/tau_g) * (u[23+k] - u[14+k]) 
    end
    
    du[22] = 0.0                             

    du[32] = (1/tau_g) * (u[34] - u[32])
    du[33] = (1/tau_g) * (u[35] - u[33])
    
    error = (Ca_tgt - Ca) 
    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH_eff, tLNS]
    for k in 0:7
        if (u[14+k] >= g_max_main[k+1] && error > 0) || (u[23+k] <= 0.0 && error < 0)
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts_main[k+1]) * error
        end
    end
    
    du[31] = 0.0
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if (u[32+k] >= g_max_cal[k+1] && error > 0) || (u[34+k] <= 0.0 && error < 0)
            du[34+k] = 0.0
        else
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end

end


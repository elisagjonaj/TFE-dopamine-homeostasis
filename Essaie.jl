"""
DA_homeo_ODE
-------------
Modèle Fyon de base avec contrôleur homéostatique intégré.

Caractéristiques :
  - Aucun blocage appliqué — simulation baseline
"""
function DA_homeo_ODE(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa = p[4:13]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    Ca = max(1e-12, u[13])

    # 3. CONDUCTANCES ET mRNA
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    
    # 5. ÉQUATION DU POTENTIEL 
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
        - gCaL * l * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        - gPace * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)

    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)

    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end
    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)

    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
end

"""
DA_homeo_ODE_blockL
--------------------

Modèle Fyon avec contrôleur homéostatique et blocage de gPace.

Caractéristiques :
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Pas de canaux N-type ni T-type
"""

function DA_homeo_ODE_blockL(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa = p[4:13]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    Ca = max(1e-12, u[13]) 

    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    gCaL_regulee = u[15] 
    gCaL_eff = (t < 30000.0) ? gPace_regulee : 1e-18
    tCaL_eff = (t < 30000.0) ? p[5] : 1e16

    # 5. ÉQUATION DU POTENTIEL
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
        - gCaL_eff * l * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        - gPace * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)

    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)

    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end
    
    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)

    ts_main = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tPace]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
end

"""
DA_homeo_ODE_blockgpace
--------------------

Modèle Fyon avec contrôleur homéostatique et blocage de gPace.

Caractéristiques :
  - Blocage instantané de gPace à t = 30000 ms (gPace_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Pas de canaux N-type ni T-type
"""

function DA_homeo_ODE_blockgpace(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa = p[4:13]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    Ca = max(1e-12, u[13]) 

    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    gPace_regulee = u[21] 
    gPace_eff = (t < 30000.0) ? gPace_regulee : 1e-18
    tPace_eff = (t < 30000.0) ? p[11] : 1e16

    # 5. ÉQUATION DU POTENTIEL
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
        - gCaL * l * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        - gPace_eff * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)

    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)

    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end
    
    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)

    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace_eff]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
end


"""
DA_homeo_2014_ODE
------------------
Modèle Yu de base avec contrôleur homéostatique intégré.

Caractéristiques :
  - Paramètres NON-physiologiques originaux du modèle Yu
    (Vhalf Na = -30.09 mV, Vhalf CaL = -45 mV)
  - Pas de canaux N-type ni T-type
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Aucun blocage appliqué — simulation baseline

"""
function DA_homeo_2014_ODE(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]   
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa = p[4:12]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances dynamiques 
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end

    
    # --- 4. DYNAMIQUE DU POTENTIEL MEMBRANAIRE ---
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gCaL * l * (V - VCa)
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

    # --- 6. DYNAMIQUE DU CALCIUM ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa)) / (F * d * 0.1)

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = (1/tau_g) * (mLCa - gLCa)

    error = (Ca_tgt - Ca) 

    ts = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa]
    for k in 0:8
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
  
end

"""
DA_homeo_2014_ODE_linf
-----------------------
Modèle Yu avec contrôleur homéostatique et shift instantané
des demi-voltages d'activation vers des valeurs physiologiques.

Caractéristiques :
  - À t = 30000 ms : Vhalf de m et l décalés vers -10 mV
  - Permet de tester si le contrôleur restaure le pacemaking
    après adoption de paramètres réalistes
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Pas de canaux N-type ni T-type
"""
function DA_homeo_2014_ODE_linf(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]   
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa = p[4:12]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end

    # --- . CALCUL DE L'ACTIVATION L-TYPE AVEC SHIFT ---
    v_half_L = (t < 30000.0) ? -45.0 : -10.0
    l_inf(V) = boltz(V, v_half_L, 7.5)
    
    v_half_N = (t < 30000.0) ? -30.09 : -10.0
    m_inf(V) = boltz(V, v_half_N, 13.2)
    
    # --- 4. DYNAMIQUE DU POTENTIEL MEMBRANAIRE ---
    du[1] = 1/C * (
        - gNa * m^3 * h * hs * (V - VNa)
        - gCaL * l * (V - VCa)
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

    # --- 6. DYNAMIQUE DU CALCIUM ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa)) / (F * d * 0.1)

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = (1/tau_g) * (mLCa - gLCa)

    error = (Ca_tgt - Ca) 

    ts = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa]
    for k in 0:8
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
  
end

"""
DA_homeo_2014_ODE_blockL
-------------------------
Modèle Yu avec contrôleur homéostatique et blocage des canaux
L-type (Cav1.3).

Caractéristiques :
  - Paramètres NON-physiologiques originaux du modèle Yu
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Pas de canaux N-type ni T-type
"""

function DA_homeo_2014_ODE_blockL(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]   
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa = p[4:12]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances 
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    

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

    # --- 6. DYNAMIQUE DU CALCIUM  ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa)) / (F * d * 0.1)

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = (1/tau_g) * (mLCa - gLCa)

    error = (Ca_tgt - Ca) 

    ts = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa]
    for k in 0:8
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
  
end

"""
DA_homeo_2014_ODE_Ntype_Ttype
------------------------------
Modèle Yu augmenté avec canaux N-type (Cav2.2) et T-type (Cav3.1),
avec contrôleur homéostatique intégré.

Caractéristiques :
  - Paramètres NON-physiologiques originaux du modèle Yu
  - Canaux N-type et T-type actifs comme sources calciques
    supplémentaires
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Aucun blocage appliqué — simulation baseline augmentée

"""
function DA_homeo_2014_ODE_Ntype_Ttype(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]    
    tau_g      = p[3]   
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa, tNCa, tTCA = p[4:14]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances 
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    
    gs_cal = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[34:35])
    
    gNtype, gTtype = gs_cal 
    mNtype, mTtype = ms_cal
    
    mN, hN, mT, hT = u[36:39]

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    
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
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa)+ gTtype * mT * hT * (V - VCa))/ (F * d * 0.1)

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = (1/tau_g) * (mLCa - gLCa)
    du[32] = (1/tau_g) * (mNtype - gNtype)
    du[33] = (1/tau_g) * (mTtype - gTtype)

    error = (Ca_tgt - Ca) 

    ts = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa]
    for k in 0:8
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0 && error <0
            du[34+k] =0.0
        else 
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end

end

"""
DA_homeo_ODE_Ntype_Ttype_gpace
--------------------------------
Modèle Fyon augmenté avec canaux N-type et T-type,
avec contrôleur homéostatique intégré.

Caractéristiques :
  - Paramètres physiologiques (Vhalf = -10 mV)
  - gPace actif, canaux N-type et T-type actifs
  - Fuite calcique (gLCa) régulée librement par le contrôleur
  - Aucun blocage appliqué — simulation baseline augmentée
"""

function DA_homeo_ODE_Ntype_Ttype_gpace(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa, tNCa, tTCA = p[4:15]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    mN, hN, mT, hT = u[38:41]
    Ca = max(1e-12, u[13]) 


    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    gs_cal = max.(0.0, u[34:35])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[36:37])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak
    gNtype, gTtype = gs_cal

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak
    mNtype, mTtype = ms_cal

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
  
    # 5. ÉQUATION DU POTENTIEL 
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
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
        - gPace * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[38] = (mN_inf(V) - mN) / tau_mN(V)
    du[39] = (hN_inf(V) - hN) / tau_hN(V)    
    du[40] = (mT_inf(V) - mT) / tau_mT(V)
    du[41] = (hT_inf(V) - hT) / tau_hT(V)


    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)
    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end

    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)
    
    du[34] = (1/tau_g) * (mNtype - gNtype)
    du[35] = (1/tau_g) * (mTtype - gTtype)

    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0.0 && error < 0
            du[36+k] = 0.0
        else
            du[36+k] = (1/ts_cal[k+1]) * error
        end
    end
           
end
"""
DA_homeo_2014_ODE_Ntype_Ttype_linfok
--------------------------------------
Modèle Yu augmenté avec canaux N-type et T-type, et blocage
des canaux L-type (Cav1.3).

Caractéristiques :
  - Paramètres NON-physiologiques originaux du modèle Yu
  - Canaux N-type et T-type actifs
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
"""
function DA_homeo_2014_ODE_Ntype_Ttype_linfok(du, u, p, t)
    # Paramètres
    Iapp       = p[1](t)
    Ca_tgt     = p[2]    
    tau_g      = p[3]    
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa, tNCa, tTCA = p[4:14]

    # --- 2. VARIABLES D'ÉTAT ---
    V, m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[1:12]
    Ca = max(1e-12, u[13]) 
    
    # Conductances 
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    
    gs_cal = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[34:35])
    
    gNtype, gTtype = gs_cal 
    mNtype, mTtype = ms_cal
    
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

    # --- 6. DYNAMIQUE DU CALCIUM  ---
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa)+ gTtype * mT * hT * (V - VCa))/ (F * d * 0.1)

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = (1/tau_g) * (mLCa - gLCa)
    du[32] = (1/tau_g) * (mNtype - gNtype)
    du[33] = (1/tau_g) * (mTtype - gTtype)

    error = (Ca_tgt - Ca) 
    ts = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tLNS, tLCa]
    for k in 0:8
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0 && error <0
            du[34+k] =0.0
        else 
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end

end

"""
DA_homeo_ODE_Ntype_Ttype_gpace_Ltype
--------------------------------------
Modèle Fyon augmenté avec canaux N-type et T-type, et blocage
des canaux L-type (Cav1.3).

Caractéristiques :
  - gPace actif, canaux N-type et T-type actifs
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
"""

function DA_homeo_ODE_Ntype_Ttype_gpace_Ltype(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa, tNCa, tTCA = p[4:15]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    mN, hN, mT, hT = u[38:41]
    Ca = max(1e-12, u[13]) 


    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    gs_cal = max.(0.0, u[34:35])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[36:37])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak
    gNtype, gTtype = gs_cal

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak
    mNtype, mTtype = ms_cal

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
   

    # 5. ÉQUATION DU POTENTIEL 
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
        - gCaL_eff * l * (V - VCa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        - gPace * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[38] = (mN_inf(V) - mN) / tau_mN(V)
    du[39] = (hN_inf(V) - hN) / tau_hN(V)    
    du[40] = (mT_inf(V) - mT) / tau_mT(V)
    du[41] = (hT_inf(V) - hT) / tau_hT(V)


    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)

    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end
   
    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)
    du[34] = (1/tau_g) * (mNtype - gNtype)
    du[35] = (1/tau_g) * (mTtype - gTtype)


    ts_main = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tPace]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0.0 && error < 0
            du[36+k] = 0.0
        else
            du[36+k] = (1/ts_cal[k+1]) * error
        end
    end
           
end
"""
test_gpace
-----------
Modèle Fyon augmenté avec canaux N-type et T-type, et blocage
de gPace — fuite calcique libre.

Caractéristiques :
  - Canaux N-type et T-type actifs
  - Blocage instantané de gPace à t = 30000 ms (gPace_eff → 1e-18)
  - Fuite calcique (gLCa) régulée librement par le contrôleur
"""
function test_gpace(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa, tNCa, tTCA = p[4:15]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    mN, hN, mT, hT = u[38:41]
    Ca = max(1e-12, u[13]) 

    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    gs_cal = max.(0.0, u[34:35])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[36:37])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak
    gNtype, gTtype = gs_cal

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak
    mNtype, mTtype = ms_cal

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    
    gPace_regulee = u[21] 
    gPace_eff = (t < 30000.0) ? gPace_regulee : 1e-18
    tPace_eff = (t < 30000.0) ? p[11] : 1e16


    # 5. ÉQUATION DU POTENTIEL
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
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
        - gPace_eff * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[38] = (mN_inf(V) - mN) / tau_mN(V)
    du[39] = (hN_inf(V) - hN) / tau_hN(V)    
    du[40] = (mT_inf(V) - mT) / tau_mT(V)
    du[41] = (hT_inf(V) - hT) / tau_hT(V)


    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)
    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end

    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = (1/tau_g) * (mLCa - gLCa)
 
    du[34] = (1/tau_g) * (mNtype - gNtype)
    du[35] = (1/tau_g) * (mTtype - gTtype)

    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace_eff]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end

    ts_leak = [tLNS, tLCa]
    for k in 0:1
        if ms_leak[k+1] <= 0.0 && error < 0
            du[32+k] = 0.0
        else
            du[32+k] = (1/ts_leak[k+1]) * error
        end
    end
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0.0 && error < 0
            du[36+k] = 0.0
        else
            du[36+k] = (1/ts_cal[k+1]) * error
        end
    end
           
end

"""
DA_homeo_ODE_Ntype_Ttype_FixedLeak_Ltype
-----------
Modèle Fyon augmenté avec canaux N-type et T-type, et blocage
de Cav1.3 — fuite calcique fixé

Caractéristiques :
  - Canaux N-type et T-type actifs
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Fuite calcique (gLCa) fixé et plus régulé par le contrôleur
"""
function DA_homeo_ODE_Ntype_Ttype_FixedLeak_Ltype(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa, tNCa, tTCA = p[4:15]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    mN, hN, mT, hT = u[38:41]
    Ca = max(1e-12, u[13]) 

    # 3. CONDUCTANCES ET mRNA
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    gs_cal  = max.(0.0, u[34:35])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])
    ms_cal  = max.(0.0, u[36:37])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak 
    gNtype, gTtype = gs_cal
    
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
   
    
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
        - gCaL_eff * l * (V - VCa)
        - gNtype * mN^2 * hN * (V - VCa)
        - gTtype * mT * hT * (V - VCa)
        - gKd * n^3 * (V - VK)
        - gKA * p_ka * (q1/2 + q2/2) * (V - VK)
        - gKERG * o * (V - VK)
        - gKSK * (V - VK) * SK_inf
        - gH * mH^2 * (V - VH)
        - gLCa * (V - VCa)
        - gLNS * (V - VLNS)
        - gPace * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_true(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[38] = (mN_inf(V) - mN) / tau_mN(V)
    du[39] = (hN_inf(V) - hN) / tau_hN(V)    
    du[40] = (mT_inf(V) - mT) / tau_mT(V)
    du[41] = (hT_inf(V) - hT) / tau_hT(V)


    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL_eff * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa))/(F * d * 0.1)
    

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)
    for k in 0:7 du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1]) end
    
    du[30] = (1/tau_g) * (ms_leak[1] - gLNS) 
    du[31] = 0.0                            

    du[34] = (1/tau_g) * (ms_cal[1] - gNtype)
    du[35] = (1/tau_g) * (ms_cal[2] - gTtype)

    ts_main = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tPace]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end  
    end
  
    if ms_leak[1] <= 0.0 && error < 0
        du[32] = 0.0
    else
        du[32] = (1/tLNS) * error
    end
                         
    du[33] = 0.0                           
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0.0 && error < 0
            du[36+k] = 0.0
        else
            du[36+k] = (1/ts_cal[k+1]) * error
        end
    end
        
end

"""
DA_homeo_2014_ODE_Ntype_Ttype_Ltypeblock_LeakC
-----------------------------------------------
Modèle Yu augmenté — version finale

Caractéristiques :
  - Paramètres NON-physiologiques originaux du modèle Yu
  - Canaux N-type et T-type actifs
  - Blocage instantané de gCaL à t = 30000 ms (gCaL_eff → 1e-18)
  - Conductance de fuite calcique (gLCa) FIXE — découplée du
    contrôleur homéostatique (du[22] = 0)
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
    gs = max.(0.0, u[14:22])
    ms = max.(0.0, u[23:31])
    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gLNS, gLCa = gs
    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mLNS, mLCa = ms
    
    gs_cal = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[34:35])
    
    gNtype, gTtype = gs_cal 
    mNtype, mTtype = ms_cal
    
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

    # --- 7. CONTRÔLEUR HOMÉOSTATIQUE ---
    du[14] = (1/tau_g) * (mNa - gNa)
    du[15] = (1/tau_g) * (mCaL - gCaL)
    du[16] = (1/tau_g) * (mKd - gKd)
    du[17] = (1/tau_g) * (mKA - gKA)
    du[18] = (1/tau_g) * (mKERG - gKERG)
    du[19] = (1/tau_g) * (mKSK - gKSK)
    du[20] = (1/tau_g) * (mH_rna - gH)
    du[21] = (1/tau_g) * (mLNS - gLNS)
    du[22] = 0.0
    du[32] = (1/tau_g) * (mNtype - gNtype)
    du[33] = (1/tau_g) * (mTtype - gTtype)


    error = (Ca_tgt - Ca) 

    ts = [t_Na, tCaL_eff, tKd, tKA, tKERG, tKSK, tH, tLNS]
    for k in 0:7
        if ms[k+1] <= 0 && error < 0
            du[23+k] = 0.0
        else
            du[23+k] = (1/ts[k+1]) * error
        end
    end
    
    du[31] = 0.0
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0 && error <0
            du[34+k] =0.0
        else 
            du[34+k] = (1/ts_cal[k+1]) * error
        end
    end
  
end

"""
DA_homeo_ODE_Ntype_Ttype_Ltype_leakC
--------------------------------------
Modèle Fyon augmenté avec canaux N-type et T-type, blocage de
gPace et conductance de fuite calcique fixe.

Caractéristiques :
  - Canaux N-type et T-type actifs
  - Blocage instantané de gPace à t = 10000 ms (gPace_eff → 1e-18)
  - Conductance de fuite calcique (gLCa) FIXE (du[31] = 0)
"""

function DA_homeo_ODE_Ntype_Ttype_Ltype_leakC(du, u, p, t)
    # 1. PARAMÈTRES
    Iapp   = p[1](t)
    Ca_tgt = p[2]
    tau_g  = p[3]
    # Constantes de temps 
    t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace, tLNS, tLCa, tNCa, tTCA = p[4:15]

    # 2. VARIABLES D'ÉTATS
    V = u[1]
    m, h, hs, l, n, p_ka, q1, q2, o, i, mH = u[2:12]
    mN, hN, mT, hT = u[38:41]
    Ca = max(1e-12, u[13]) 


    # 3. CONDUCTANCES ET mRNA 
    gs_main = max.(0.0, u[14:21])
    gs_leak = max.(0.0, u[30:31])
    gs_cal = max.(0.0, u[34:35])
    
    ms_main = max.(0.0, u[22:29])
    ms_leak = max.(0.0, u[32:33])
    ms_cal = max.(0.0, u[36:37])

    gNa, gCaL, gKd, gKA, gKERG, gKSK, gH, gPace = gs_main
    gLNS, gLCa = gs_leak
    gNtype, gTtype = gs_cal

    mNa, mCaL, mKd, mKA, mKERG, mKSK, mH_rna, mPace = ms_main
    mLNS, mLCa = ms_leak
    mNtype, mTtype = ms_cal

    # --- 3. DÉPENDANCE AU CALCIUM (SK & Pompe) ---
    SK_inf = 0.0
    ICap = 0.0
    if Ca > 0
        SK_inf = 1 / (1 + (0.00019 / Ca)^4)
        ICap = ICapmax / (1 + (0.0005 / Ca))
    end
    
    l_inf_local(V)  = boltz(V, -10.0, 7.5)
    gPace_regulee = u[21] 
    gPace_eff = (t < 10000.0) ? gPace_regulee : 1e-18
    tPace_eff = (t < 10000.0) ? p[11] : 1e16


    # 5. ÉQUATION DU POTENTIEL 
    du[1] = 1/C * (
        - gNa * m^3 * h * (V - VNa)
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
        - gPace_eff * mPacemaker_inf(V) * (V - EPacemaker)
        + 100 * Iapp / (pi * d * L)
    )

    # 6. DYNAMIQUE DES PORTES 
    du[2]  = (1 / tau_m(V))  * (m_inf_true(V) - m)
    du[3]  = (1 / tau_h(V))  * (h_inf(V)      - h)
    du[4]  = (1 / tau_hs(V)) * (hs_inf(V)     - hs)
    du[5]  = (1 / tau_l(V))  * (l_inf_local(V) - l)
    du[6]  = (1 / tau_n(V))  * (n_inf(V)      - n)
    du[7]  = (1 / tau_p(V))  * (p_inf(V)      - p_ka)
    du[8]  = (1 / tau_q1(V)) * (q1_inf(V)     - q1)
    du[9]  = (1 / tau_q2(V)) * (q2_inf(V)     - q2)
    du[10] = alphao(V) * (1 - o - i) + betai(V) * i - o * (alphai(V) + betao(V))
    du[11] = alphai(V) * o - betai(V) * i
    du[12] = (1 / tau_mH(V)) * (mH_inf(V) - mH)
    du[38] = (mN_inf(V) - mN) / tau_mN(V)
    du[39] = (hN_inf(V) - hN) / tau_hN(V)    
    du[40] = (mT_inf(V) - mT) / tau_mT(V)
    du[41] = (hT_inf(V) - hT) / tau_hT(V)


    # 7. DYNAMIQUE DU CALCIUM 
    du[13] = -2 * fCa * (gLCa * (V - VCa) + ICap + gCaL * l * (V - VCa) + gNtype * mN^2 * hN * (V - VCa) + gTtype * mT * hT * (V - VCa)) / (F * d * 0.1)

    # 8. CONTRÔLEUR HOMÉOSTATIQUE 
    error = (Ca_tgt - Ca)
    
    for k in 0:7
        du[14+k] = (1/tau_g) * (ms_main[k+1] - gs_main[k+1])
    end
    
    du[30] = (1/tau_g) * (mLNS - gLNS)
    du[31] = 0.0
  
    du[34] = (1/tau_g) * (mNtype - gNtype)
    du[35] = (1/tau_g) * (mTtype - gTtype)

    ts_main = [t_Na, tCaL, tKd, tKA, tKERG, tKSK, tH, tPace_eff]
    for k in 0:7
        if ms_main[k+1] <= 0.0 && error < 0
            du[22+k] = 0.0
        else
            du[22+k] = (1/ts_main[k+1]) * error
        end
    end
    
    if ms_leak[1] <= 0.0 && error < 0
        du[32] = 0.0
    else
        du[32] = (1/tLNS) * error
    end
    
    du[33] = 0.0
    
    ts_cal = [tNCa, tTCA]
    for k in 0:1
        if ms_cal[k+1] <= 0.0 && error < 0
            du[36+k] = 0.0
        else
            du[36+k] = (1/ts_cal[k+1]) * error
        end
    end
           
end
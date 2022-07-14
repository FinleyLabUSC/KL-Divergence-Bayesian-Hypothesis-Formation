
def simulatejakstat(t,y,params):
    # note-input parameters should be length 61, input species should be length 56
    PRL = y[0]
    RJ = y[1]
    S5Ac = y[2]
    S5Bc = y[3]
    SHP2 = y[4]
    PPX = y[5]
    PPN = y[6]
    PRLRJ = y[7]
    PRLRJ2 = y[8]
    PRLRJ2a = y[9]
    PRLRJ2aS5Ac = y[10]
    PRLRJ2aS5Bc = y[11]
    pS5Ac = y[12]
    pS5Bc = y[13]
    pS5AcpS5Ac = y[14]
    pS5AcpS5Bc = y[15]
    pS5BcpS5Bc = y[16]
    PRLRJ2aSHP2 = y[17]
    PPXpS5Ac = y[18]
    PPXpS5Bc = y[19]
    PPXpS5AcpS5Ac = y[20]
    PPXpS5BcpS5Bc = y[21]
    PPXpS5AcpS5Bc = y[22]
    pS5AcS5Ac = y[23]
    pS5BcS5Bc = y[24]
    S5AcpS5Bc = y[25]
    pS5AcS5Bc = y[26]
    pS5AnpS5An = y[27]
    pS5AnpS5Bn = y[28]
    pS5BnpS5Bn = y[29]
    pS5An = y[30]
    pS5Bn = y[31]
    PPNpS5An = y[32]
    PPNpS5Bn = y[33]
    S5An = y[34]
    S5Bn = y[35]
    PPNpS5AnpS5An = y[36]
    PPNpS5AnpS5Bn = y[37]
    PPNpS5BnpS5Bn = y[38]
    pS5AnS5An = y[39]
    pS5BnS5Bn = y[40]
    S5AnpS5Bn = y[41]
    pS5AnS5Bn = y[42]
    mRNAn = y[43]
    mRNAc = y[44]
    SOCS1 = y[45]
    SOCS1PRLRJ2a = y[46]
    PRLRJ2aS5AcSHP2 = y[47]
    PRLRJ2aS5BcSHP2 = y[48]
    SOCS1PRLRJ2aSHP2 = y[49]
    mRn = y[50]
    mRc = y[51]
    Rc = y[52]
    mBCLn = y[53]
    mBCLc = y[54]
    BCL = y[55]


    k1 = params[0]
    k2 = params[1]
    k_2 = params[2]
    k3 = params[3]
    k_3 = params[4]
    k4 = params[5]
    k5 = params[6]
    k_5 = params[7]
    k6 = params[8]
    kdeg = params[9]
    deg_ratio = params[10]
    k8A = params[11]
    k_8A = params[12]
    k8B = params[13]
    k_8B = params[14]
    k8AB = params[15]
    k_8AB = params[16]
    k9 = params[17]
    k_9 = params[18]
    k10 = params[19]
    k11 = params[20]
    k_11 = params[21]
    k12 = params[22]
    k13 = params[23]
    k_13 = params[24]
    k14A = params[25]
    k14B = params[26]
    k14AB = params[27]
    k15 = params[28]
    k_15 = params[29]
    k16 = params[30]
    k17inA = params[31]
    k17outA = params[32]
    k17inB = params[33]
    k17outB = params[34]
    k18a = params[35]
    k18b = params[36]
    k19 = params[37]
    k20 = params[38]
    k21 = params[39]
    k_21 = params[40]
    k22 = params[41]
    k23 = params[42]
    k24 = params[43]
    k25a = params[44]
    k25b = params[45]
    k26 = params[46]
    k27 = params[47]
    k28 = params[48]
    k29 = params[49]
    Vratio = params[50]
    ncratioA = params[51]
    ncratioB = params[52]
    totalSTAT = params[53]
    k30a = params[54]
    k30b = params[55]
    k31 = params[56]
    k32 = params[57]
    k33 = params[58]
    k34 = params[59]
    k35 = params[60]


    dPRLdt = (k_2 * PRLRJ - k2 * PRL * RJ) * 1.39E-4 - (1 * 0.693 / (6 * 3600) * PRL)
    '#Adjusted for volume of cell to volume of cytoplasm ratio, 1st order deg with 6 hour half life'
    dRJdt = k1 + k_2 * PRLRJ - k2 * PRL * RJ - k3 * PRLRJ * RJ + k_3 * PRLRJ2 - kdeg * RJ + k29 * Rc
    dS5Acdt = k_5 * PRLRJ2aS5Ac - k5 * S5Ac * PRLRJ2a + k12 * PPXpS5Ac - k13 * S5Ac * pS5Ac + k_13 * pS5AcS5Ac\
              - k13 * S5Ac * pS5Bc + k_13 * S5AcpS5Bc + k10 * PRLRJ2aS5AcSHP2 + kdeg * deg_ratio\
              * (PRLRJ2aS5Ac + PRLRJ2aS5AcSHP2) - k17inA * S5Ac + k17outA * S5An * Vratio
    dS5Bcdt = k_5 * PRLRJ2aS5Bc - k5 * S5Bc * PRLRJ2a + k12 * PPXpS5Bc - k13 * S5Bc * pS5Bc + k_13\
              * pS5BcS5Bc - k13 * S5Bc * pS5Ac + k_13 * pS5AcS5Bc + k10 * PRLRJ2aS5BcSHP2 + kdeg * deg_ratio\
              * (PRLRJ2aS5Bc + PRLRJ2aS5BcSHP2) - k17inB * S5Bc + k17outB * S5Bn * Vratio
    dSHP2dt = k_9 * PRLRJ2aSHP2 - k9 * PRLRJ2a * SHP2 + k10 * PRLRJ2aSHP2 - k9 * SHP2 * SOCS1PRLRJ2a + k_9\
              * SOCS1PRLRJ2aSHP2 + k10 * SOCS1PRLRJ2aSHP2 - k9 * SHP2 * PRLRJ2aS5Ac + k_9 * PRLRJ2aS5AcSHP2\
              - k9 * SHP2 * PRLRJ2aS5Bc + k_9 * PRLRJ2aS5BcSHP2 + k10 * PRLRJ2aS5AcSHP2 + k10 * PRLRJ2aS5BcSHP2\
              + kdeg * deg_ratio * (PRLRJ2aSHP2 + PRLRJ2aS5AcSHP2 + PRLRJ2aS5BcSHP2 + SOCS1PRLRJ2aSHP2)\
              + k24 * SOCS1PRLRJ2aSHP2
    dPPXdt = k_11 * PPXpS5Ac - k11 * PPX * pS5Ac + k12 * PPXpS5Ac + k_11 * PPXpS5Bc - k11 * PPX * pS5Bc\
                + k12 * PPXpS5Bc - k11 * PPX * pS5AcpS5Ac + k_11 * PPXpS5AcpS5Ac + k12 * PPXpS5AcpS5Ac\
                - k11 * PPX * pS5BcpS5Bc + k_11 * PPXpS5BcpS5Bc + k12 * PPXpS5BcpS5Bc - k11 * PPX * pS5AcpS5Bc\
                + k_11 * PPXpS5AcpS5Bc + k12 * PPXpS5AcpS5Bc
    dPPNdt = k_15 * PPNpS5An - k15 * PPN * pS5An + k16 * PPNpS5An + k_15 * PPNpS5Bn - k15 * PPN * pS5Bn\
                + k16 * PPNpS5Bn - k15 * PPN * pS5AnpS5An + k_15 * PPNpS5AnpS5An + k16 * PPNpS5AnpS5An\
                - k15 * PPN * pS5BnpS5Bn + k_15 * PPNpS5BnpS5Bn + k16 * PPNpS5BnpS5Bn - k15 * PPN * pS5AnpS5Bn\
                + k_15 * PPNpS5AnpS5Bn + k16 * PPNpS5AnpS5Bn
    dPRLRJdt = k2 * PRL * RJ - k_2 * PRLRJ - k3 * PRLRJ * RJ + k_3 * PRLRJ2 - kdeg * deg_ratio * PRLRJ
    dPRLRJ2dt = k3 * PRLRJ * RJ - k_3 * PRLRJ2 - k4 * PRLRJ2 + k10 * PRLRJ2aSHP2 + k10 * SOCS1PRLRJ2aSHP2\
                + k10 * PRLRJ2aS5AcSHP2 + k10 * PRLRJ2aS5BcSHP2 - kdeg * deg_ratio * PRLRJ2
    dPRLRJ2adt = k4 * PRLRJ2 - k5 * PRLRJ2a * S5Ac + k_5 * PRLRJ2aS5Ac + k6 * PRLRJ2aS5Ac - k5 * PRLRJ2a * S5Bc\
                 + k_5 * PRLRJ2aS5Bc + k6 * PRLRJ2aS5Bc - k9 * PRLRJ2a * SHP2 + k_9 * PRLRJ2aSHP2\
                 - k21 * PRLRJ2a * SOCS1 + k_21 * SOCS1PRLRJ2a + k23 * SOCS1PRLRJ2a - kdeg * deg_ratio * PRLRJ2a
    dPRLRJ2aS5Acdt = k5 * PRLRJ2a * S5Ac - k_5 * PRLRJ2aS5Ac - k6 * PRLRJ2aS5Ac - k9 * SHP2 * PRLRJ2aS5Ac\
                 + k_9 * PRLRJ2aS5AcSHP2 - kdeg * deg_ratio * PRLRJ2aS5Ac
    dPRLRJ2aS5Bcdt = k5 * PRLRJ2a * S5Bc - k_5 * PRLRJ2aS5Bc - k6 * PRLRJ2aS5Bc - k9 * SHP2 * PRLRJ2aS5Bc\
                 + k_9 * PRLRJ2aS5BcSHP2 - kdeg * deg_ratio * PRLRJ2aS5Bc
    dpS5Acdt = k6 * PRLRJ2aS5Ac - 2 * k8A * pS5Ac * pS5Ac + 2 * k_8A * pS5AcpS5Ac - k8AB * pS5Ac * pS5Bc\
                 + k_8AB * pS5AcpS5Bc - k11 * PPX * pS5Ac + k_11 * PPXpS5Ac - k13 * pS5Ac * S5Ac\
                 + k_13 * pS5AcS5Ac - k13 * pS5Ac * S5Bc + k_13 * pS5AcS5Bc
    dpS5Bcdt = k6 * PRLRJ2aS5Bc - 2 * k8B * pS5Bc * pS5Bc + 2 * k_8B * pS5BcpS5Bc - k8AB * pS5Ac * pS5Bc\
                 + k_8AB * pS5AcpS5Bc - k11 * PPX * pS5Bc + k_11 * PPXpS5Bc - k13 * pS5Bc * S5Bc + k_13 * pS5BcS5Bc\
                 - k13 * pS5Bc * S5Ac + k_13 * S5AcpS5Bc
    dpS5AcpS5Acdt = k8A * pS5Ac * pS5Ac - k_8A * pS5AcpS5Ac - k11 * PPX * pS5AcpS5Ac + k_11 * PPXpS5AcpS5Ac\
                 - k14A * pS5AcpS5Ac
    dpS5AcpS5Bcdt = k8AB * pS5Ac * pS5Bc - k_8AB * pS5AcpS5Bc - k11 * PPX * pS5AcpS5Bc + k_11 * PPXpS5AcpS5Bc\
                 - k14AB * pS5AcpS5Bc
    dpS5BcpS5Bcdt = k8B * pS5Bc * pS5Bc - k_8B * pS5BcpS5Bc - k11 * PPX * pS5BcpS5Bc + k_11 * PPXpS5BcpS5Bc\
                 - k14B * pS5BcpS5Bc
    dPRLRJ2aSHP2dt = k9 * PRLRJ2a * SHP2 - k_9 * PRLRJ2aSHP2 - k10 * PRLRJ2aSHP2 + k23 * SOCS1PRLRJ2aSHP2\
                 - kdeg * deg_ratio * PRLRJ2aSHP2
    dPPXpS5Acdt = k11 * PPX * pS5Ac - k_11 * PPXpS5Ac - k12 * PPXpS5Ac
    dPPXpS5Bcdt = k11 * PPX * pS5Bc - k_11 * PPXpS5Bc - k12 * PPXpS5Bc
    dPPXpS5AcpS5Acdt = k11 * PPX * pS5AcpS5Ac - k_11 * PPXpS5AcpS5Ac - k12 * PPXpS5AcpS5Ac
    dPPXpS5BcpS5Bcdt = k11 * PPX * pS5BcpS5Bc - k_11 * PPXpS5BcpS5Bc - k12 * PPXpS5BcpS5Bc
    dPPXpS5AcpS5Bcdt = k11 * PPX * pS5AcpS5Bc - k_11 * PPXpS5AcpS5Bc - k12 * PPXpS5AcpS5Bc
    dpS5AcS5Acdt = k13 * pS5Ac * S5Ac - k_13 * pS5AcS5Ac + k12 * PPXpS5AcpS5Ac
    dpS5BcS5Bcdt = k13 * pS5Bc * S5Bc - k_13 * pS5BcS5Bc + k12 * PPXpS5BcpS5Bc
    dS5AcpS5Bcdt = k13 * S5Ac * pS5Bc - k_13 * S5AcpS5Bc + 0.5 * k12 * PPXpS5AcpS5Bc
    dpS5AcS5Bcdt = k13 * pS5Ac * S5Bc - k_13 * pS5AcS5Bc + 0.5 * k12 * PPXpS5AcpS5Bc
    dpS5AnpS5Andt = k14A * pS5AcpS5Ac/Vratio + k8A * pS5An * pS5An - k_8A * pS5AnpS5An - k15 * PPN * pS5AnpS5An\
                 + k_15 * PPNpS5AnpS5An
    dpS5AnpS5Bndt = k14AB * pS5AcpS5Bc/Vratio + k8AB * pS5An * pS5Bn - k_8AB * pS5AnpS5Bn - k15 * PPN * pS5AnpS5Bn\
                 + k_15 * PPNpS5AnpS5Bn
    dpS5BnpS5Bndt = k14B * pS5BcpS5Bc/Vratio + k8B * pS5Bn * pS5Bn - k_8B * pS5BnpS5Bn - k15 * PPN * pS5BnpS5Bn\
                 + k_15 * PPNpS5BnpS5Bn
    dpS5Andt = 2 * k_8A * pS5AnpS5An - 2 * k8A * pS5An * pS5An + k_8AB * pS5AnpS5Bn - k8AB * pS5An * pS5Bn\
                 - k15 * PPN * pS5An + k_15 * PPNpS5An - k13 * pS5An * S5An + k_13 * pS5AnS5An - k13 * pS5An\
                 * S5Bn + k_13 * pS5AnS5Bn
    dpS5Bndt = 2 * k_8B * pS5BnpS5Bn - 2 * k8B * pS5Bn * pS5Bn + k_8AB * pS5AnpS5Bn - k8AB * pS5An * pS5Bn\
                 - k15 * PPN * pS5Bn + k_15 * PPNpS5Bn - k13 * pS5Bn * S5Bn + k_13 * pS5BnS5Bn - k13 * pS5Bn\
                 * S5An + k_13 * S5AnpS5Bn
    dPPNpS5Andt = k15 * PPN * pS5An - k_15 * PPNpS5An - k16 * PPNpS5An
    dPPNpS5Bndt = k15 * PPN * pS5Bn - k_15 * PPNpS5Bn - k16 * PPNpS5Bn
    dS5Andt = k16 * PPNpS5An - k13 * pS5An * S5An + k_13 * pS5AnS5An - k13 * pS5Bn * S5An + k_13 * S5AnpS5Bn\
                 + k17inA * S5Ac/Vratio - k17outA * S5An
    dS5Bndt = k16 * PPNpS5Bn - k13 * pS5Bn * S5Bn + k_13 * pS5BnS5Bn - k13 * pS5An * S5Bn + k_13 * pS5AnS5Bn\
                 + k17inB * S5Bc/Vratio - k17outB * S5Bn
    dPPNpS5AnpS5Andt = k15 * PPN * pS5AnpS5An - k_15 * PPNpS5AnpS5An - k16 * PPNpS5AnpS5An
    dPPNpS5AnpS5Bndt = k15 * PPN * pS5AnpS5Bn - k_15 * PPNpS5AnpS5Bn - k16 * PPNpS5AnpS5Bn
    dPPNpS5BnpS5Bndt = k15 * PPN * pS5BnpS5Bn - k_15 * PPNpS5BnpS5Bn - k16 * PPNpS5BnpS5Bn
    dpS5AnS5Andt = k16 * PPNpS5AnpS5An + k13 * pS5An * S5An - k_13 * pS5AnS5An
    dpS5BnS5Bndt = k16 * PPNpS5BnpS5Bn + k13 * pS5Bn * S5Bn - k_13 * pS5BnS5Bn
    dS5AnpS5Bndt = 0.5 * k16 * PPNpS5AnpS5Bn + k13 * S5An * pS5Bn - k_13 * S5AnpS5Bn
    dpS5AnS5Bndt = 0.5 * k16 * PPNpS5AnpS5Bn + k13 * pS5An * S5Bn - k_13 * pS5AnS5Bn
    dmRNAndt = k18a * (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn) / (k18b + (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn))\
                 - k19 * mRNAn
    dmRNAcdt = k19 * mRNAn* Vratio - k22 * mRNAc
    dSOCS1dt = k20 * mRNAc - k21 * SOCS1 * PRLRJ2a + k_21 * SOCS1PRLRJ2a - k23 * SOCS1 + k10 * SOCS1PRLRJ2aSHP2\
                 + kdeg * deg_ratio * (SOCS1PRLRJ2a + SOCS1PRLRJ2aSHP2)
    dSOCS1PRLRJ2adt = k21 * SOCS1 * PRLRJ2a - k_21 * SOCS1PRLRJ2a - k9 * SHP2 * SOCS1PRLRJ2a + k_9 * SOCS1PRLRJ2aSHP2\
                 - k23 * SOCS1PRLRJ2a - kdeg * deg_ratio * SOCS1PRLRJ2a - k24 * SOCS1PRLRJ2a
    dPRLRJ2aS5AcSHP2dt = k9 * SHP2 * PRLRJ2aS5Ac - k_9 * PRLRJ2aS5AcSHP2 - k10 * PRLRJ2aS5AcSHP2 - kdeg * deg_ratio\
                 * PRLRJ2aS5AcSHP2
    dPRLRJ2aS5BcSHP2dt = k9 * SHP2 * PRLRJ2aS5Bc - k_9 * PRLRJ2aS5BcSHP2 - k10 * PRLRJ2aS5BcSHP2 - kdeg * deg_ratio\
                 * PRLRJ2aS5BcSHP2
    dSOCS1PRLRJ2aSHP2dt = k9 * SHP2 * SOCS1PRLRJ2a - k_9 * SOCS1PRLRJ2aSHP2 - k10 * SOCS1PRLRJ2aSHP2 - k23 * SOCS1PRLRJ2aSHP2\
                 - kdeg * deg_ratio * SOCS1PRLRJ2aSHP2 - k24 * SOCS1PRLRJ2aSHP2
    dmRndt = k25a * (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn) / (k25b + (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn))\
                 - k26 * mRn
    dmRcdt = k26 * mRn * Vratio - k27 * mRc
    dRcdt = k28 * mRc - k29 * Rc
    dmBCLndt = k30a * (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn) / (k30b + (pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn))\
                 - k31 * mBCLn
    dmBCLcdt = k31 * mBCLn * Vratio - k32 * mBCLc
    dBCLdt = k33 * mBCLc - k34 * BCL + k35
    return [dPRLdt,dRJdt,dS5Acdt,dS5Bcdt,dSHP2dt,dPPXdt,dPPNdt,dPRLRJdt,dPRLRJ2dt,dPRLRJ2adt,dPRLRJ2aS5Acdt,
            dPRLRJ2aS5Bcdt,dpS5Acdt,dpS5Bcdt,dpS5AcpS5Acdt,dpS5AcpS5Bcdt,dpS5BcpS5Bcdt,dPRLRJ2aSHP2dt, dPPXpS5Acdt,
            dPPXpS5Bcdt,dPPXpS5AcpS5Acdt,dPPXpS5BcpS5Bcdt,dPPXpS5AcpS5Bcdt,dpS5AcS5Acdt,dpS5BcS5Bcdt,dS5AcpS5Bcdt,
            dpS5AcS5Bcdt,dpS5AnpS5Andt,dpS5AnpS5Bndt,dpS5BnpS5Bndt,dpS5Andt,dpS5Bndt,dPPNpS5Andt,dPPNpS5Bndt,dS5Andt,
            dS5Bndt,dPPNpS5AnpS5Andt,dPPNpS5AnpS5Bndt,dPPNpS5BnpS5Bndt,dpS5AnS5Andt,dpS5BnS5Bndt,dS5AnpS5Bndt,
            dpS5AnS5Bndt,dmRNAndt,dmRNAcdt,dSOCS1dt,dSOCS1PRLRJ2adt,dPRLRJ2aS5AcSHP2dt,dPRLRJ2aS5BcSHP2dt,
            dSOCS1PRLRJ2aSHP2dt,dmRndt,dmRcdt,dRcdt,dmBCLndt,dmBCLcdt,dBCLdt]


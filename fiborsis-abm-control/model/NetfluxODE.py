import numpy as np
from scipy.integrate import odeint

def ODEfunc(t, y, tau, ymax, w, n, EC50):
    AngII = 0
    AT1R = 1
    AGT = 2
    ACE = 3
    NOX = 4
    ROS = 5
    ET1 = 6
    ETAR = 7
    DAG = 8
    PKC = 9
    TRPC = 10
    NE = 11
    BAR = 12
    Forskolin = 13
    AC = 14
    cAMP = 15
    PKA = 16
    CREB = 17
    CBP = 18
    TGFB = 19
    TGFB1R = 20
    smad3 = 21
    smad7 = 22
    latentTGFB = 23
    BAMBI = 24
    PDGF = 25
    PDGFR = 26
    NP = 27
    NPRA = 28
    cGMP = 29
    PKG = 30
    mechanical = 31
    B1int = 32
    Rho = 33
    ROCK = 34
    Ca = 35
    calcineurin = 36
    NFAT = 37
    IL6 = 38
    gp130 = 39
    STAT = 40
    IL1 = 41
    IL1RI = 42
    TNFa = 43
    TNFaR = 44
    NFKB = 45
    PI3K = 46
    Akt = 47
    p38 = 48
    TRAF = 49
    ASK1 = 50
    MKK3 = 51
    PP1 = 52
    JNK = 53
    abl = 54
    Rac1 = 55
    MEKK1 = 56
    MKK4 = 57
    ERK = 58
    Ras = 59
    Raf = 60
    MEK1 = 61
    FAK = 62
    epac = 63
    Factin = 64
    FA = 65
    migration = 66
    cmyc = 67
    CTGF = 68
    proliferation = 69
    SRF = 70
    EDAFN = 71
    aSMA = 72
    AP1 = 73
    TIMP1 = 74
    TIMP2 = 75
    PAI1 = 76
    proMMP14 = 77
    proMMP1 = 78
    proMMP2 = 79
    proMMP9 = 80
    MMP1 = 81
    MMP2 = 82
    MMP9 = 83
    MMP14 = 84
    fibronectin = 85
    periostin = 86
    CImRNA = 87
    CIIImRNA = 88
    CI = 89
    CIII = 90

    dydt = np.zeros(91)
    dydt[AngII] = (OR(w[0], AND(w[13], [act(y[AGT], w[13], n[13], EC50[13]), act(y[ACE], w[13], n[13], EC50[13])])) *
                   ymax[AngII] - y[AngII]) / tau[AngII]
    dydt[AT1R] = (act(y[AngII], w[18], n[18], EC50[18]) * ymax[AT1R] - y[AT1R]) / tau[AT1R]
    dydt[AGT] = (AND(w[30], [inhib(y[AT1R], w[30], n[30], EC50[30]), act(y[p38], w[30], n[30], EC50[30]),
                             inhib(y[JNK], w[30], n[30], EC50[30])]) * ymax[AGT] - y[AGT]) / tau[AGT]
    dydt[ACE] = (act(y[TGFB1R], w[53], n[53], EC50[53]) * ymax[ACE] - y[ACE]) / tau[ACE]
    dydt[NOX] = (OR(act(y[AT1R], w[19], n[19], EC50[19]), act(y[TGFB1R], w[116], n[116], EC50[116])) * ymax[NOX] - y[
        NOX]) / tau[NOX]
    dydt[ROS] = (OR(act(y[NOX], w[20], n[20], EC50[20]), act(y[ETAR], w[40], n[40], EC50[40])) * ymax[ROS] - y[ROS]) / \
                tau[ROS]
    dydt[ET1] = (OR(w[8], act(y[AP1], w[17], n[17], EC50[17])) * ymax[ET1] - y[ET1]) / tau[ET1]
    dydt[ETAR] = (act(y[ET1], w[65], n[65], EC50[65]) * ymax[ETAR] - y[ETAR]) / tau[ETAR]
    dydt[DAG] = (OR(act(y[ETAR], w[136], n[136], EC50[136]), act(y[AT1R], w[137], n[137], EC50[137])) * ymax[DAG] - y[
        DAG]) / tau[DAG]
    dydt[PKC] = (act(y[DAG], w[141], n[141], EC50[141]) * ymax[PKC] - y[PKC]) / tau[PKC]
    dydt[TRPC] = (act(y[DAG], w[138], n[138], EC50[138]) * ymax[TRPC] - y[TRPC]) / tau[TRPC]
    dydt[NE] = (w[6] * ymax[NE] - y[NE]) / tau[NE]
    dydt[BAR] = (act(y[NE], w[64], n[64], EC50[64]) * ymax[BAR] - y[BAR]) / tau[BAR]
    dydt[Forskolin] = (w[10] * ymax[Forskolin] - y[Forskolin]) / tau[Forskolin]
    dydt[AC] = (OR(act(y[BAR], w[77], n[77], EC50[77]),
                   OR(AND(w[78], [act(y[AT1R], w[78], n[78], EC50[78]), act(y[BAR], w[78], n[78], EC50[78])]),
                      act(y[Forskolin], w[121], n[121], EC50[121]))) * ymax[AC] - y[AC]) / tau[AC]
    dydt[cAMP] = (act(y[AC], w[79], n[79], EC50[79]) * ymax[cAMP] - y[cAMP]) / tau[cAMP]
    dydt[PKA] = (act(y[cAMP], w[47], n[47], EC50[47]) * ymax[PKA] - y[PKA]) / tau[PKA]
    dydt[CREB] = (act(y[PKA], w[62], n[62], EC50[62]) * ymax[CREB] - y[CREB]) / tau[CREB]
    dydt[CBP] = (OR(inhib(y[smad3], w[49], n[49], EC50[49]), inhib(y[CREB], w[50], n[50], EC50[50])) * ymax[CBP] - y[
        CBP]) / tau[CBP]
    dydt[TGFB] = (OR(w[1],
                     OR(AND(w[11], [act(y[latentTGFB], w[11], n[11], EC50[11]), act(y[MMP9], w[11], n[11], EC50[11])]),
                        AND(w[12],
                            [act(y[latentTGFB], w[12], n[12], EC50[12]), act(y[MMP2], w[12], n[12], EC50[12])]))) *
                  ymax[TGFB] - y[TGFB]) / tau[TGFB]
    dydt[TGFB1R] = (AND(w[60], [act(y[TGFB], w[60], n[60], EC50[60]), inhib(y[BAMBI], w[60], n[60], EC50[60])]) * ymax[
        TGFB1R] - y[TGFB1R]) / tau[TGFB1R]
    dydt[smad3] = (AND(w[31], [act(y[TGFB1R], w[31], n[31], EC50[31]), inhib(y[smad7], w[31], n[31], EC50[31]),
                               inhib(y[PKG], w[31], n[31], EC50[31])]) * ymax[smad3] - y[smad3]) / tau[smad3]
    dydt[smad7] = (act(y[STAT], w[122], n[122], EC50[122]) * ymax[smad7] - y[smad7]) / tau[smad7]
    dydt[latentTGFB] = (act(y[AP1], w[81], n[81], EC50[81]) * ymax[latentTGFB] - y[latentTGFB]) / tau[latentTGFB]
    dydt[BAMBI] = (AND(w[120], [act(y[TGFB], w[120], n[120], EC50[120]), act(y[IL1RI], w[120], n[120], EC50[120])]) *
                   ymax[BAMBI] - y[BAMBI]) / tau[BAMBI]
    dydt[PDGF] = (w[7] * ymax[PDGF] - y[PDGF]) / tau[PDGF]
    dydt[PDGFR] = (act(y[PDGF], w[76], n[76], EC50[76]) * ymax[PDGFR] - y[PDGFR]) / tau[PDGFR]
    dydt[NP] = (w[9] * ymax[NP] - y[NP]) / tau[NP]
    dydt[NPRA] = (act(y[NP], w[87], n[87], EC50[87]) * ymax[NPRA] - y[NPRA]) / tau[NPRA]
    dydt[cGMP] = (act(y[NPRA], w[88], n[88], EC50[88]) * ymax[cGMP] - y[cGMP]) / tau[cGMP]
    dydt[PKG] = (act(y[cGMP], w[89], n[89], EC50[89]) * ymax[PKG] - y[PKG]) / tau[PKG]
    dydt[mechanical] = (w[2] * ymax[mechanical] - y[mechanical]) / tau[mechanical]
    dydt[B1int] = (OR(AND(w[46], [act(y[PKC], w[46], n[46], EC50[46]), act(y[mechanical], w[46], n[46], EC50[46])]),
                      act(y[mechanical], w[51], n[51], EC50[51])) * ymax[B1int] - y[B1int]) / tau[B1int]
    dydt[Rho] = (act(y[B1int], w[27], n[27], EC50[27]) * ymax[Rho] - y[Rho]) / tau[Rho]
    dydt[ROCK] = (act(y[Rho], w[83], n[83], EC50[83]) * ymax[ROCK] - y[ROCK]) / tau[ROCK]
    dydt[Ca] = (act(y[TRPC], w[139], n[139], EC50[139]) * ymax[Ca] - y[Ca]) / tau[Ca]
    dydt[calcineurin] = (act(y[Ca], w[140], n[140], EC50[140]) * ymax[calcineurin] - y[calcineurin]) / tau[calcineurin]
    dydt[NFAT] = (act(y[calcineurin], w[131], n[131], EC50[131]) * ymax[NFAT] - y[NFAT]) / tau[NFAT]
    dydt[IL6] = (OR(w[3], OR(AND(w[14], [act(y[CREB], w[14], n[14], EC50[14]), act(y[CBP], w[14], n[14], EC50[14])]),
                             OR(act(y[NFKB], w[15], n[15], EC50[15]), act(y[AP1], w[16], n[16], EC50[16])))) * ymax[
                     IL6] - y[IL6]) / tau[IL6]
    dydt[gp130] = (act(y[IL6], w[21], n[21], EC50[21]) * ymax[gp130] - y[gp130]) / tau[gp130]
    dydt[STAT] = (act(y[gp130], w[28], n[28], EC50[28]) * ymax[STAT] - y[STAT]) / tau[STAT]
    dydt[IL1] = (w[4] * ymax[IL1] - y[IL1]) / tau[IL1]
    dydt[IL1RI] = (act(y[IL1], w[67], n[67], EC50[67]) * ymax[IL1RI] - y[IL1RI]) / tau[IL1RI]
    dydt[TNFa] = (w[5] * ymax[TNFa] - y[TNFa]) / tau[TNFa]
    dydt[TNFaR] = (act(y[TNFa], w[86], n[86], EC50[86]) * ymax[TNFaR] - y[TNFaR]) / tau[TNFaR]
    dydt[NFKB] = (OR(act(y[IL1RI], w[25], n[25], EC50[25]), OR(act(y[ERK], w[37], n[37], EC50[37]),
                                                               OR(act(y[p38], w[38], n[38], EC50[38]),
                                                                  act(y[Akt], w[117], n[117], EC50[117])))) * ymax[
                      NFKB] - y[NFKB]) / tau[NFKB]
    dydt[PI3K] = (OR(act(y[TNFaR], w[29], n[29], EC50[29]), OR(act(y[TGFB1R], w[113], n[113], EC50[113]),
                                                               OR(act(y[PDGFR], w[114], n[114], EC50[114]),
                                                                  act(y[FAK], w[115], n[115], EC50[115])))) * ymax[
                      PI3K] - y[PI3K]) / tau[PI3K]
    dydt[Akt] = (act(y[PI3K], w[112], n[112], EC50[112]) * ymax[Akt] - y[Akt]) / tau[Akt]
    dydt[p38] = (OR(act(y[ROS], w[23], n[23], EC50[23]), OR(act(y[MKK3], w[94], n[94], EC50[94]),
                                                            OR(act(y[Ras], w[111], n[111], EC50[111]), AND(w[124], [
                                                                act(y[Rho], w[124], n[124], EC50[124]),
                                                                inhib(y[Rac1], w[124], n[124], EC50[124])])))) * ymax[
                     p38] - y[p38]) / tau[p38]
    dydt[TRAF] = (OR(act(y[TGFB1R], w[95], n[95], EC50[95]), act(y[TNFaR], w[103], n[103], EC50[103])) * ymax[TRAF] - y[
        TRAF]) / tau[TRAF]
    dydt[ASK1] = (OR(act(y[TRAF], w[104], n[104], EC50[104]), act(y[IL1RI], w[107], n[107], EC50[107])) * ymax[ASK1] -
                  y[ASK1]) / tau[ASK1]
    dydt[MKK3] = (act(y[ASK1], w[105], n[105], EC50[105]) * ymax[MKK3] - y[MKK3]) / tau[MKK3]
    dydt[PP1] = (act(y[p38], w[93], n[93], EC50[93]) * ymax[PP1] - y[PP1]) / tau[PP1]
    dydt[JNK] = (OR(act(y[ROS], w[24], n[24], EC50[24]),
                    OR(AND(w[98], [inhib(y[NFKB], w[98], n[98], EC50[98]), act(y[MKK4], w[98], n[98], EC50[98])]),
                       AND(w[125],
                           [inhib(y[Rho], w[125], n[125], EC50[125]), act(y[MKK4], w[125], n[125], EC50[125])]))) *
                 ymax[JNK] - y[JNK]) / tau[JNK]
    dydt[abl] = (act(y[PDGFR], w[99], n[99], EC50[99]) * ymax[abl] - y[abl]) / tau[abl]
    dydt[Rac1] = (OR(act(y[B1int], w[26], n[26], EC50[26]), act(y[abl], w[100], n[100], EC50[100])) * ymax[Rac1] - y[
        Rac1]) / tau[Rac1]
    dydt[MEKK1] = (OR(act(y[FAK], w[80], n[80], EC50[80]), act(y[Rac1], w[96], n[96], EC50[96])) * ymax[MEKK1] - y[
        MEKK1]) / tau[MEKK1]
    dydt[MKK4] = (OR(act(y[MEKK1], w[97], n[97], EC50[97]), act(y[ASK1], w[106], n[106], EC50[106])) * ymax[MKK4] - y[
        MKK4]) / tau[MKK4]
    dydt[ERK] = (OR(act(y[ROS], w[22], n[22], EC50[22]),
                    AND(w[92], [inhib(y[PP1], w[92], n[92], EC50[92]), act(y[MEK1], w[92], n[92], EC50[92])])) * ymax[
                     ERK] - y[ERK]) / tau[ERK]
    dydt[Ras] = (act(y[AT1R], w[132], n[132], EC50[132]) * ymax[Ras] - y[Ras]) / tau[Ras]
    dydt[Raf] = (act(y[Ras], w[90], n[90], EC50[90]) * ymax[Raf] - y[Raf]) / tau[Raf]
    dydt[MEK1] = (AND(w[91], [inhib(y[ERK], w[91], n[91], EC50[91]), act(y[Raf], w[91], n[91], EC50[91])]) * ymax[
        MEK1] - y[MEK1]) / tau[MEK1]
    dydt[FAK] = (act(y[ROCK], w[133], n[133], EC50[133]) * ymax[FAK] - y[FAK]) / tau[FAK]
    dydt[epac] = (act(y[cAMP], w[82], n[82], EC50[82]) * ymax[epac] - y[epac]) / tau[epac]
    dydt[Factin] = (act(y[ROCK], w[126], n[126], EC50[126]) * ymax[Factin] - y[Factin]) / tau[Factin]
    dydt[FA] = (AND(w[128], [act(y[B1int], w[128], n[128], EC50[128]), act(y[Factin], w[128], n[128], EC50[128])]) *
                ymax[FA] - y[FA]) / tau[FA]
    dydt[migration] = (OR(act(y[MMP9], w[84], n[84], EC50[84]), OR(act(y[MMP2], w[85], n[85], EC50[85]), AND(w[109], [
        inhib(y[PKA], w[109], n[109], EC50[109]), act(y[epac], w[109], n[109], EC50[109])]))) * ymax[migration] - y[
                           migration]) / tau[migration]
    dydt[cmyc] = (act(y[JNK], w[101], n[101], EC50[101]) * ymax[cmyc] - y[cmyc]) / tau[cmyc]
    dydt[CTGF] = (AND(w[32], [act(y[CBP], w[32], n[32], EC50[32]), act(y[smad3], w[32], n[32], EC50[32]),
                              act(y[ERK], w[32], n[32], EC50[32])]) * ymax[CTGF] - y[CTGF]) / tau[CTGF]
    dydt[proliferation] = (OR(act(y[AP1], w[61], n[61], EC50[61]), OR(act(y[CREB], w[63], n[63], EC50[63]),
                                                                      OR(act(y[CTGF], w[66], n[66], EC50[66]),
                                                                         OR(act(y[PKC], w[68], n[68], EC50[68]),
                                                                            act(y[cmyc], w[102], n[102],
                                                                                EC50[102]))))) * ymax[proliferation] -
                           y[proliferation]) / tau[proliferation]
    dydt[SRF] = (act(y[Factin], w[127], n[127], EC50[127]) * ymax[SRF] - y[SRF]) / tau[SRF]
    dydt[EDAFN] = (act(y[NFAT], w[52], n[52], EC50[52]) * ymax[EDAFN] - y[EDAFN]) / tau[EDAFN]
    dydt[aSMA] = (OR(AND(w[130], [act(y[CBP], w[130], n[130], EC50[130]), act(y[smad3], w[130], n[130], EC50[130]),
                                  act(y[SRF], w[130], n[130], EC50[130])]),
                     OR(AND(w[134], [act(y[CBP], w[134], n[134], EC50[134]), act(y[smad3], w[134], n[134], EC50[134])]),
                        act(y[SRF], w[135], n[135], EC50[135]))) * ymax[aSMA] - y[aSMA]) / tau[aSMA]
    dydt[AP1] = (OR(act(y[ERK], w[41], n[41], EC50[41]), act(y[JNK], w[119], n[119], EC50[119])) * ymax[AP1] - y[AP1]) / \
                tau[AP1]
    dydt[TIMP1] = (act(y[AP1], w[44], n[44], EC50[44]) * ymax[TIMP1] - y[TIMP1]) / tau[TIMP1]
    dydt[TIMP2] = (act(y[AP1], w[45], n[45], EC50[45]) * ymax[TIMP2] - y[TIMP2]) / tau[TIMP2]
    dydt[PAI1] = (act(y[smad3], w[108], n[108], EC50[108]) * ymax[PAI1] - y[PAI1]) / tau[PAI1]
    dydt[proMMP14] = (OR(act(y[AP1], w[75], n[75], EC50[75]), act(y[NFKB], w[110], n[110], EC50[110])) * ymax[
        proMMP14] - y[proMMP14]) / tau[proMMP14]
    dydt[proMMP1] = (AND(w[39], [inhib(y[smad3], w[39], n[39], EC50[39]), act(y[NFKB], w[39], n[39], EC50[39]),
                                 act(y[AP1], w[39], n[39], EC50[39])]) * ymax[proMMP1] - y[proMMP1]) / tau[proMMP1]
    dydt[proMMP2] = (OR(act(y[STAT], w[33], n[33], EC50[33]), act(y[AP1], w[42], n[42], EC50[42])) * ymax[proMMP2] - y[
        proMMP2]) / tau[proMMP2]
    dydt[proMMP9] = (OR(act(y[STAT], w[34], n[34], EC50[34]),
                        AND(w[43], [act(y[NFKB], w[43], n[43], EC50[43]), act(y[AP1], w[43], n[43], EC50[43])])) * ymax[
                         proMMP9] - y[proMMP9]) / tau[proMMP9]
    dydt[MMP1] = (AND(w[56], [inhib(y[TIMP1], w[56], n[56], EC50[56]), act(y[proMMP1], w[56], n[56], EC50[56])]) * ymax[
        MMP1] - y[MMP1]) / tau[MMP1]
    dydt[MMP2] = (OR(AND(w[58], [inhib(y[TIMP1], w[58], n[58], EC50[58]), act(y[proMMP2], w[58], n[58], EC50[58]),
                                 act(y[MMP14], w[58], n[58], EC50[58])]),
                     AND(w[59], [act(y[proMMP2], w[59], n[59], EC50[59]), act(y[MMP14], w[59], n[59], EC50[59])])) *
                  ymax[MMP2] - y[MMP2]) / tau[MMP2]
    dydt[MMP9] = (OR(AND(w[55], [inhib(y[TIMP1], w[55], n[55], EC50[55]), act(y[proMMP9], w[55], n[55], EC50[55])]),
                     AND(w[57], [inhib(y[TIMP2], w[57], n[57], EC50[57]), act(y[proMMP9], w[57], n[57], EC50[57])])) *
                  ymax[MMP9] - y[MMP9]) / tau[MMP9]
    dydt[MMP14] = (act(y[proMMP14], w[54], n[54], EC50[54]) * ymax[MMP14] - y[MMP14]) / tau[MMP14]
    dydt[fibronectin] = (OR(AND(w[48], [act(y[CBP], w[48], n[48], EC50[48]), act(y[smad3], w[48], n[48], EC50[48])]),
                            act(y[NFKB], w[118], n[118], EC50[118])) * ymax[fibronectin] - y[fibronectin]) / tau[
                            fibronectin]
    dydt[periostin] = (OR(AND(w[35], [act(y[CBP], w[35], n[35], EC50[35]), act(y[smad3], w[35], n[35], EC50[35])]),
                          AND(w[36], [act(y[CREB], w[36], n[36], EC50[36]), act(y[CBP], w[36], n[36], EC50[36])])) *
                       ymax[periostin] - y[periostin]) / tau[periostin]
    dydt[CImRNA] = (OR(AND(w[69], [act(y[CBP], w[69], n[69], EC50[69]), act(y[smad3], w[69], n[69], EC50[69]),
                                   inhib(y[epac], w[69], n[69], EC50[69])]), act(y[SRF], w[123], n[123], EC50[123])) *
                    ymax[CImRNA] - y[CImRNA]) / tau[CImRNA]
    dydt[CIIImRNA] = (OR(AND(w[70], [act(y[CBP], w[70], n[70], EC50[70]), act(y[smad3], w[70], n[70], EC50[70]),
                                     inhib(y[epac], w[70], n[70], EC50[70])]), act(y[SRF], w[129], n[129], EC50[129])) *
                      ymax[CIIImRNA] - y[CIIImRNA]) / tau[CIIImRNA]
    dydt[CI] = (OR(AND(w[71], [inhib(y[MMP1], w[71], n[71], EC50[71]), act(y[CImRNA], w[71], n[71], EC50[71])]),
                   AND(w[73], [inhib(y[MMP2], w[73], n[73], EC50[73]), act(y[CImRNA], w[73], n[73], EC50[73])])) * ymax[
                    CI] - y[CI]) / tau[CI]
    dydt[CIII] = (OR(AND(w[72], [inhib(y[MMP1], w[72], n[72], EC50[72]), act(y[CIIImRNA], w[72], n[72], EC50[72])]),
                     AND(w[74], [inhib(y[MMP2], w[74], n[74], EC50[74]), act(y[CIIImRNA], w[74], n[74], EC50[74])])) *
                  ymax[CIII] - y[CIII]) / tau[CIII]

    return dydt


# utility functions

def act(x, w, n, EC50):
    # hill activation function with parameters w (weight), n (Hill coeff), EC50
    beta = ((EC50 ** n) - 1) / (2 * EC50 ** n - 1)
    K = (beta - 1) ** (1 / n)
    fact = w * (beta * np.sign(x) * np.abs(x) ** n) / (K ** n + np.abs(x) ** n)
    if fact > w:
        fact = w
    return fact


def inhib(x, w, n, EC50):
    # inverse hill function with parameters w (weight), n (Hill coeff), EC50
    finhib = w - act(x, w, n, EC50)
    return finhib


def OR(x, y):
    # OR logic gate
    z = x + y - x * y
    return z


def AND(w, reactList):
    # AND logic gate, multiplying all of the reactants together
    if w == 0:
        z = 0
    else:
        p = np.array(reactList).prod()
        z = p / w ** (len(reactList) - 2)
    return z


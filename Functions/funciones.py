# Codigo para calcular la velocidad de corrosion en pozos petroleros de acuerdo con el modelo NORSOK M-506 2017, y el indice de saturacion basado en el modelo de Oddo Tomson
#
# VARIABLES QUE INGRESAN AL MODELO
#
# Deberiamos extrar esta informacion de los archivos Excel del cliente o hacer formularios para que ingrese la informacion que falte
#
# pres: total absolute pressure in well head [psia]
# pres1: total absolute pressure in well bottom [psia]
# temp: temperature in well head [F]
# temp1: temperature in well bottom [F]
# Na: sodium concentration in water [mg/l]
# K: potasium concentration in water [mg/l]
# Mg: magnesium concentration in water [mg/l]
# Ca: calcium concentration in water [mg/l]
# Sr: strontium concentration in water [mg/l]
# Ba: barium concentration in water [mg/l]
# Cl: chloride concentration in water [mg/l]
# SO4: sulfate concentration in water [mg/l]
# HCO3: total alkalinity of water, as mg/l of HCO3
# co2fraction: mol fraction of CO2 in gas
# cac: carboxylic acids concentration as acetate in water [mg/l]
# bopd: barrels of oil per day [BOPD]
# bwpd: barrels of water per day [BWPD]
# mscf: volumetric flowrate of gas in thousand of standard cubic feet per day [MSCFD].
# diameter: pipe internal diameter [in]
# h1: well depth [ft]
# IC_efficiency: %, corrosion inhibitor efficiency

# Paquetes necesarios

# In[1]:


import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


# Funcion para calcular ionic strength

# In[2]:


def ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3):
    cna = (Na / 22990) * (1) ** 2
    ck = (K / 39100) * (1) ** 2
    cmg = (Mg / 24310) * (2) ** 2
    cca = (Ca / 40080) * (2) ** 2
    csr = (Sr / 87620) * (2) ** 2
    cba = (Ba / 137340) * (2) ** 2
    ccl = (Cl / 35450) * (-1) ** 2
    cso4 = (SO4 / 96060) * (-2) ** 2
    chco3 = (HCO3 / 61018) * (-1) ** 2
    ionic = 0.5 * (cna + ck + cmg + cca + csr + cba + ccl + cso4 + chco3)
    return ionic


# Funcion para calcular el coeficiente de fugacidad del CO2 segun Oddo-Tomson

# In[3]:


def fugacityCO2(temp, pres):
    fugacity = math.exp(pres * (2.84e-4 - (0.255 / (temp + 460))))
    return fugacity


# Funcion para calcular la fraccion del CO2

# In[4]:


def ytgas(co2fraction, bopd, bwpd, mscf, temp, pres):
    mmscf = mscf / 1000
    fg = fugacityCO2(temp, pres)
    denominador = 1 + ((pres * fg * (5 * bwpd + 10 * bopd) * 10e-6) / (mmscf * (temp + 460)))
    yg = co2fraction / denominador
    return yg


# Funcion para calcular la alcalinidad (HCO3) segun Odo-Tomson

# In[5]:


def hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf):
    # cac: molar concentration of acetic acid ion
    cac = cac / 59040
    # hco3: molar concentration of bicarbonate ion
    HCO3 = HCO3 / 61018
    fg = fugacityCO2(temp, pres)
    yg = ytgas(co2fraction, bopd, bwpd, mscf, temp, pres)
    pco2 = pres * fg * yg
    ionics = ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
    pk = 3.79 + 6.80e-3 * temp - 13.2e-6 * temp ** 2 - 3.66e-5 * pres - 0.097 * ionics ** (1 / 2) + 0.221 * ionics
    kc = 10 ** (-pk)
    # kac = 10 ** -(4.81 - 1.49e-3 * temp + 1.095e-5 * temp**2 - 1.16e-5 * pres - 0.893 * ionics**(1/2) + 0.437 * ionics)
    a = 1
    b = kc * pco2 + cac - HCO3
    c = -kc * pco2 * HCO3
    numerador = -b + ((b ** 2) - 4 * a * c) ** (1 / 2)
    hco3 = numerador / 2
    return hco3


# Funcion para calcular el pH del agua

# In[6]:


def ph(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf):
    alkalinity = hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
    si = ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
    mco2 = ytgas(co2fraction, bopd, bwpd, mscf, temp, pres)
    fugacity = fugacityCO2(temp, pres)
    k = 8.60 + 5.31e-3 * temp - 2.253e-6 * temp ** 2 - 2.237e-5 * pres - 0.990 * si ** (1 / 2) + 0.658 * si
    phcal = math.log10(alkalinity / (pres * fugacity * mco2)) + k
    return phcal


# Funcion para calcular el indice de saturacion de la calcita

# In[7]:


def indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf):
    calcium = Ca / 40080
    alkalinity = hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
    ionic = ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
    k = 5.85 + 15.19e-3 * temp - 1.64e-6 * temp ** 2 - 5.27e-5 * pres - 3.334 * ionic ** (1 / 2) + 1.431 * ionic
    fugacity = fugacityCO2(temp, pres)
    mco2 = ytgas(co2fraction, bopd, bwpd, mscf, temp, pres)
    log_result = math.log10((calcium * (alkalinity) ** 2) / (pres * fugacity * mco2))
    indice = log_result + k
    return indice


# Funcion para calcular la cantidad de escala que se forma (PTB). Basado en el paper: "Scale prediction, measurement and control in a subsea satellite field" (Wright, 1994)

# In[8]:


def ptb(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf):
    calcium = Ca / 40080
    alkalinity = hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
    si = indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
    G = calcium + alkalinity
    X = calcium - alkalinity
    K = calcium * alkalinity * (10 ** (-si))
    ptb_c = 17500 * (G - (X ** 2 + (4 * K)) ** (0.5))
    if ptb_c > 0:
        ptbc1 = ptb_c
    else:
        ptbc1 = 0
    return ptbc1


# Funcion para calcular el factor pH

# In[9]:


def fpH_FixT(tempe, iph):
    tempo = 7
    if tempe == 5.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 2.0676 - 0.2309 * iph
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 4.342 - (1.051 * iph) + (0.0708 * iph ** 2)
    if tempe == 15.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 2.0676 - (0.2309 * iph)
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 4.986 - (1.191 * iph) + (0.0708 * iph ** 2)
    if tempe == 20.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 2.0676 - (0.2309 * iph)
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 5.1885 - (1.2353 * iph) + (0.0708 * iph ** 2)
    if tempe == 40.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 2.0676 - (0.2309 * iph)
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 5.1885 - (1.2353 * iph) + (0.0708 * iph ** 2)
    if tempe == 60.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 1.836 - (0.1818 * iph)
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 15.444 - (6.1291 * iph) + (0.8204 * iph ** 2) - (0.0371 * iph ** 3)
    if tempe == 80.0:
        if (iph >= 3.5) and (iph < 4.6):
            tempo = 2.6727 - (0.3636 * iph)
        if (iph >= 4.6) and (iph <= 6.5):
            tempo = 331.68 * math.exp(-1.2618 * iph)
    if tempe == 90.0:
        if (iph >= 3.5) and (iph < 4.57):
            tempo = 3.1355 - (0.4673 * iph)
        if (iph >= 4.57) and (iph < 5.62):
            tempo = 21254 * math.exp(-2.1811 * iph)
        if (iph >= 5.62) and (iph <= 6.5):
            tempo = 0.4014 - (0.0538 * iph)
    if tempe == 120.0:
        if (iph >= 3.5) and (iph < 4.3):
            tempo = 1.5375 - (0.125 * iph)
        if (iph >= 4.3) and (iph < 5.0):
            tempo = 5.9757 - 1.157 * iph
        if (iph >= 5.0) and (iph <= 6.5):
            tempo = 0.546125 - (0.071225 * iph)
    if tempe == 150.0:
        if (iph >= 3.5) and (iph < 3.8):
            tempo = 1
        if (iph >= 3.8) and (iph < 5.0):
            tempo = 17.634 - (7.0945 * iph) + (0.715 * iph ** 2)
        if (iph >= 5.0) and (iph <= 6.5):
            tempo = 0.037

    return tempo


# Funcion para interpolar el factor pH

# In[10]:


def fpH_Cal(tempe, iph):
    TempRange = [5.0, 15.0, 20.0, 40.0, 60.0, 80.0, 90.0, 120.0, 150.0]

    loc = 0

    for i, temp_i in enumerate(TempRange):
        if temp_i > tempe:
            loc = i
            break
    TempLower = TempRange[loc - 1]
    TempUpper = TempRange[loc]

    fpHLower = fpH_FixT(TempLower, iph)
    fpHUpper = fpH_FixT(TempUpper, iph)

    tempo = (fpHUpper - fpHLower) / (TempUpper - TempLower)

    tempo = fpHLower + (tempe - TempLower) * tempo

    return tempo


# Funcion para calcular la fugacidad del CO2 segun Norsok

# In[11]:


def FugacityofCO2(co2fraction, pres, temp, bopd, bwpd, mscf):
    # co2 fraction: mole fraction of CO2 in the gas phase (<1)

    # Conversion pressure to bar
    pres_bar = pres * 0.06894

    # Conversion temperature to Kelvin
    temp_K = (temp - 32) * (5 / 9) + 273.15

    # Calculation of fugacity coefficient
    A = 10 ** (pres_bar * (0.0031 - 1.4 / temp_K))

    # Calculation of CO2 partial pressure
    co2fraction = ytgas(co2fraction, bopd, bwpd, mscf, temp, pres)

    co2pressure = co2fraction * pres_bar

    # Calculation of CO2 fugacity
    tempo = A * co2pressure

    return tempo


# Funcion para calcular el esfuerzo cortante en la tuberia

# In[12]:


def Shearstress(mscf, bopd, bwpd, diameter, temp, pres):
    # mscf: gas flowrate in thousand of standard cubic feet per day [MSCFD]
    # bopd: barrels of oil per day [BOPD]
    # bwpd: barrels of water per day [BWPD]
    # diameter: pipe internal diameter [in]
    # temp: temperature [F]
    # pres: pressure well head [psi]

    # Converting pipe diameter to m
    diameter = diameter * 0.0254

    # Converting temp to C
    temp = (temp - 32) * (5 / 9)

    # Converting pressure to bar
    pres = pres * 0.0689476

    # Converting volumetric flowrate of liquid to m3/s
    vol_l = (bopd + bwpd) * 1.84E-6

    # Converting volumetric flowrate of gas to m3/s
    vol_g = (mscf * 28.3168) / (24 * 3600)

    # Calculation of liquid and gas superficial velocity (m/s)
    v_sl = vol_l / ((3.1416 * diameter ** 2) / 4)

    v_sg = (vol_g / ((3.1416 * diameter ** 2) / 4)) * 0.9 * ((temp + 273.15) / (288.7))

    # Calculation of superficial velocity of mixture
    v_m = v_sg + v_sl

    # Calculation of liquid fraction (holdup)
    holdup = vol_l / (vol_l + vol_g)

    # density of water [kg/m3]. Default value in Norsok: 1024 kg/m3
    den_w = 998.2 / (1 + (0.0002 * (temp - 20)))

    # density of oil [kg/m3]. Default value in Norsok: 850 kg/m3
    den_o = 850 - ((1.825 - (0.001315 * 850)) * (temp - 20))

    # Calculation density of liquid
    water_cut = bwpd / (bwpd + bopd)
    density_l = (water_cut * den_w) + (den_o * (1 - water_cut))

    # Calculation of density of gas. Default in Norsok: 10 kg/m3
    # density_g=28.97*((pres*0.75)/(0.9*8.314*(temp+273.25)*1000))
    density_g = (2.7 * 14.5 * 16.018 * pres * 0.5537) / (0.9 * (460 + (1.8 * temp + 32)))

    # Calculation of density of mixture
    density_mix = density_l * holdup + density_g * (1 - holdup)

    # Viscosity of gas conversion from cp to kg/m.s. I use the default value of Norsok (given in cp) for viscosity of gas.
    vis_g = 0.03 / 1000

    # Viscosity of oil from https://www.petroskills.com/blog/entry/crude-oil-and-changing-temperature#.X8i1prlxc2w for a crude of 25 API
    # Viscosity equation from "Estimating the Viscosity of Crude Oil System", SPE5434
    sg = 0.9042 - (5.93E-4) * (temp - 15)
    api = (141.5 / sg) - 131.5
    y = 10 ** (3.0324 - 0.02023 * api)
    x = y * (1.8 * temp + 32) ** -1.163
    vis_o = ((10 ** x) - 1) / 1000

    # Viscosity of water
    vis_w = ((temp + 273.15) - 225.4) ** -1.637

    # Calculation viscosity of liquid
    if water_cut < 0.5:
        bo1 = water_cut / (1.187 * (1 - (1 / 7.06) ** 0.4))
        vis_l = vis_o * ((1 + ((water_cut / bo1) / (1.187 - (water_cut / bo1)))) ** 2.5)
    else:
        bo2 = (1 - water_cut) / (1.187 * (1 - ((vis_w / vis_o) / 7.06) ** 0.4))
        vis_l = vis_w * (1 + (((1 - water_cut) / bo2) / (1.187 - ((1 - water_cut) / bo2)))) ** 2.5

    # Calculation viscosity of mixture
    vis_m = vis_l * holdup + vis_g * (1 - holdup)

    # Calculation friction factor. I use the default value of roughness [m] of Norsok
    roughness = 50E-6
    friction = 0.001375 * (
            1 + ((20000 * roughness / diameter) + (10 ** 6 * (vis_m / (v_m * diameter * density_mix)))) ** 0.33)

    # Calculation of shear stress
    tempo = 0.5 * density_mix * friction * v_m ** 2
    return tempo


# Funcion para calcular la constante Kt utilizada en el modelo NORSOK

# In[13]:


def Kt(temp):
    # Converting temp to C
    temp = (temp - 32) * (5 / 9)

    temp_table = (5, 15, 20, 40, 60, 80, 90, 120, 150)
    value_table = (0.42, 1.59, 4.762, 8.927, 10.695, 9.949, 6.250, 7.770, 5.203)

    for i, temp_i in enumerate(temp_table):

        if temp <= temp_i:
            temp_lower = temp_table[i - 1]
            temp_upper = temp_table[i]
            value_lower = value_table[i - 1]
            value_upper = value_table[i]
            break

    return round(value_upper - (temp_upper - temp) * (value_upper - value_lower) / (temp_upper - temp_lower), 3)


# Funcion para calcular la velocidad de corrosion

# In[14]:


def calcNorsok(temp, pres, bopd, bwpd, mscf, co2fraction, HCO3, Cl, Na, K, Mg, Ca, Sr, Ba, SO4, cac, diameter,
               IC_efficiency, situacion):
    phcal = ph(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
    fphcal = fpH_Cal((temp - 32) * (5 / 9), phcal)
    fy = FugacityofCO2(co2fraction, pres, temp, bopd, bwpd, mscf)
    kt = Kt(temp)
    shs = Shearstress(mscf, bopd, bwpd, diameter, temp, pres)

    # Calculo de la velocidad de corrosion segun Norsok, en mm/year
    nk = kt * fy ** 0.62 * (shs / 19) ** (0.146 + 0.0324 * math.log10(fy)) * fphcal

    # Transformar a mpy y tomar en cuenta la eficiencia del inhibidor de corrosion
    corr_ic = nk * 39.4 * (100 - IC_efficiency) / 100

    # Asignar niveles de riesgo segun norma NACE
    if corr_ic < 1:
        corr_risk = 'Bajo'
    if corr_ic >= 1 and corr_ic < 5:
        corr_risk = 'Moderado'
    if corr_ic >= 5 and corr_ic < 10:
        corr_risk = "Alto"
    if corr_ic >= 10:
        corr_risk = 'Muy alto'

    return nk, corr_ic, corr_risk


# Funcion para calcular el indice de saturacion

# In[15]:


def calCalcite(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf, situacion):
    calcite_si = indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)

    solid = ptb(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)

    # Asignar niveles de riesgo
    if calcite_si <= 0.5:
        scale_risk = 'Bajo'
    if calcite_si > 0.5 and calcite_si <= 1.5:
        scale_risk = 'Moderado'
    if calcite_si > 1.5 and calcite_si <= 2.5:
        scale_risk = 'Alto'
    if calcite_si > 2.5:
        scale_risk = 'Muy alto'

    return calcite_si, solid, scale_risk


# Funcion para graficar el perfil de la velocidad de corrosion

# In[16]:


def grahpNorskok(temp, temp1, pres, pres1, bopd, bwpd, mscf, co2fraction, HCO3, Cl, Na, K, Mg, Ca, Sr, Ba, SO4, cac,
                 diameter, IC_efficiency, h1, well_name):
    temp_array = np.linspace(temp, temp1, 100)
    press_array = np.linspace(pres, pres1, 100)
    depth_array = np.linspace(0, h1, 100)
    aux = 0
    kt_df = []
    yt_df = []
    fy_df = []
    shearstress_df = []
    ph_df = []
    fph_df = []
    nk_df = []
    corr_profile_risk = []

    for i in temp_array:
        auxkt = Kt(i)
        kt_df.append(auxkt)
        auxyt = ytgas(co2fraction, bopd, bwpd, mscf, i, press_array[aux])
        yt_df.append(auxyt)
        auxfy = FugacityofCO2(co2fraction, press_array[aux], i, bopd, bwpd, mscf)
        fy_df.append(auxfy)
        auxss = Shearstress(mscf, bopd, bwpd, diameter, i, press_array[aux])
        shearstress_df.append(auxss)
        auxph = ph(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
        ph_df.append(auxph)
        auxfph = fpH_Cal((i - 32) * (5 / 9), auxph)
        fph_df.append(auxfph)
        auxnk = auxkt * auxfy ** 0.62 * (auxss / 19) ** (0.146 + 0.0324 * math.log10(auxfy)) * auxfph * 39.4 * (
                100 - IC_efficiency) / 100
        nk_df.append(auxnk)

        # Asignar niveles de riesgo segun norma NACE
        if auxnk < 1:
            auxcorr = 'Bajo'
        if auxnk >= 1 and auxnk < 5:
            auxcorr = 'Moderado'
        if auxnk >= 5 and auxnk < 10:
            auxcorr = "Alto"
        if auxnk >= 10:
            auxcorr = 'Muy alto'

        corr_profile_risk.append(auxcorr)

        aux = aux + 1

    return temp_array, press_array, depth_array, fy_df, ph_df, nk_df, corr_profile_risk


# Funcion para graficar el perfil del indice de saturacion

# In[17]:


def graphCalcite(temp, temp1, pres, pres1, h1, bopd, bwpd, mscf, co2fraction, HCO3, Cl, Na, K, Mg, Ca, Sr, Ba, SO4, cac,
                 well_name):
    temp_array = np.linspace(temp, temp1, 100)
    press_array = np.linspace(pres, pres1, 100)
    depth_array = np.linspace(0, h1, 100)
    aux = 0
    alkalinity = []
    fy = []
    yt = []
    ph1 = []
    calcite = []
    ptb1 = []
    scale_profile_risk = []

    for i in temp_array:
        auxak = hco(press_array[aux], i, Na, K, Mg, Ca, Sr,
                    Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
        alkalinity.append(auxak)
        auxfy = fugacityCO2(i, press_array[aux])
        fy.append(auxfy)
        auxyt = ytgas(co2fraction, bopd, bwpd, mscf, i, press_array[aux])
        yt.append(auxyt)
        auxph = ph(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
        ph1.append(auxph)
        auxcalcite = indiceSaturacion(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd,
                                      bwpd, mscf)
        calcite.append(auxcalcite)
        auxptb = ptb(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
        ptb1.append(auxptb)

        # Asignar niveles de riesgo
        if auxcalcite <= 0.5:
            auxscale = 'Bajo'
        if auxcalcite > 0.5 and auxcalcite <= 1.5:
            auxscale = 'Moderado'
        if auxcalcite > 1.5 and auxcalcite <= 2.5:
            auxscale = 'Alto'
        if auxcalcite > 2.5:
            auxscale = 'Muy alto'

        scale_profile_risk.append(auxscale)

        aux = aux + 1

    return temp_array, press_array, depth_array, fy, ph1, calcite, ptb1, scale_profile_risk

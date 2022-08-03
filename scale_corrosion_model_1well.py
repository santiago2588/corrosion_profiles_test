#!/usr/bin/env python
# coding: utf-8

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
    cna = (Na / 22990) * (1)  ** 2
    ck = (K / 39100)  * (1) ** 2
    cmg = (Mg / 24310) * (2) ** 2
    cca = (Ca / 40080) * (2) ** 2
    csr = (Sr / 87620) * (2) ** 2
    cba = (Ba / 137340) * (2) ** 2
    ccl = (Cl / 35450) * (-1) ** 2
    cso4 = (SO4 / 96060) * (-2) ** 2
    chco3 = (HCO3 / 61018) * (-1) **2
    ionic = 0.5 * (cna + ck + cmg + cca + csr + cba + ccl + cso4 + chco3)
    return ionic


# Funcion para calcular el coeficiente de fugacidad del CO2 segun Oddo-Tomson

# In[3]:


def fugacityCO2(temp, pres):
  fugacity = math.exp(pres * (2.84e-4 - (0.255 / (temp + 460))))
  return fugacity


# Funcion para calcular la fraccion del CO2

# In[4]:


def ytgas(co2fraction,bopd,bwpd,mscf,temp,pres):
    mmscf=mscf/1000
    fg = fugacityCO2(temp,pres)
    denominador = 1 + ((pres * fg * (5 * bwpd + 10 * bopd) * 10e-6)/ (mmscf * (temp + 460))) 
    yg = co2fraction / denominador
    return yg


# Funcion para calcular la alcalinidad (HCO3) segun Odo-Tomson

# In[5]:


def hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf):
    #cac: molar concentration of acetic acid ion
    cac=cac/59040
    #hco3: molar concentration of bicarbonate ion
    HCO3=HCO3/61018
    fg=fugacityCO2(temp,pres)
    yg=ytgas(co2fraction,bopd,bwpd,mscf,temp,pres)
    pco2 = pres * fg * yg
    ionics=ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
    pk=3.79 + 6.80e-3 * temp - 13.2e-6 * temp**2 - 3.66e-5 * pres - 0.097 * ionics**(1/2) + 0.221 * ionics
    kc = 10 **(-pk)
    #kac = 10 ** -(4.81 - 1.49e-3 * temp + 1.095e-5 * temp**2 - 1.16e-5 * pres - 0.893 * ionics**(1/2) + 0.437 * ionics)
    a = 1
    b = kc * pco2 + cac - HCO3
    c = -kc * pco2 * HCO3
    numerador = -b + ((b**2) - 4 * a* c)**(1/2)
    hco3 = numerador / 2 
    return hco3


# Funcion para calcular el pH del agua

# In[6]:


def ph(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf):
  alkalinity=hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
  si=ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
  mco2=ytgas(co2fraction,bopd,bwpd,mscf,temp,pres)
  fugacity=fugacityCO2(temp,pres)
  k = 8.60 + 5.31e-3 * temp - 2.253e-6 * temp**2 - 2.237e-5 * pres - 0.990 * si**(1/2) + 0.658 * si
  phcal = math.log10(alkalinity/(pres * fugacity * mco2)) + k
  return phcal


# Funcion para calcular el indice de saturacion de la calcita

# In[7]:


def indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf):
  calcium=Ca/40080
  alkalinity=hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
  ionic=ionicS(Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3)
  k = 5.85 + 15.19e-3 * temp - 1.64e-6 * temp**2 - 5.27e-5 * pres - 3.334 * ionic**(1/2) + 1.431 * ionic
  fugacity=fugacityCO2(temp,pres)
  mco2=ytgas(co2fraction,bopd,bwpd,mscf,temp,pres)
  log_result = math.log10((calcium * (alkalinity)**2) / (pres * fugacity * mco2))
  indice = log_result + k
  return indice


# Funcion para calcular la cantidad de escala que se forma (PTB). Basado en el paper: "Scale prediction, measurement and control in a subsea satellite field" (Wright, 1994)

# In[8]:


def ptb(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf):
    calcium=Ca/40080
    alkalinity=hco(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
    si=indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
    G=calcium+alkalinity
    X=calcium-alkalinity
    K=calcium*alkalinity*(10**(-si))
    ptb_c=17500*(G-(X**2+(4*K))**(0.5))
    if ptb_c > 0: 
        ptbc1=ptb_c
    else:
        ptbc1=0
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
    pres_bar = pres*0.06894
    
    #Conversion temperature to Kelvin
    temp_K = (temp-32)*(5/9)+273.15
    
    #Calculation of fugacity coefficient
    A = 10 ** (pres_bar * (0.0031 - 1.4 / temp_K))
    
    #Calculation of CO2 partial pressure
    co2fraction= ytgas(co2fraction,bopd,bwpd,mscf,temp,pres)
         
    co2pressure = co2fraction * pres_bar
    
    #Calculation of CO2 fugacity
    tempo = A * co2pressure
    
    return tempo


# Funcion para calcular el esfuerzo cortante en la tuberia

# In[12]:


def Shearstress(mscf,bopd,bwpd,diameter,temp,pres):
    #mscf: gas flowrate in thousand of standard cubic feet per day [MSCFD]
    #bopd: barrels of oil per day [BOPD]
    #bwpd: barrels of water per day [BWPD]
    #diameter: pipe internal diameter [in]
    #temp: temperature [F]
    #pres: pressure well head [psi]
    
    #Converting pipe diameter to m
    diameter=diameter*0.0254
    
    #Converting temp to C
    temp=(temp-32)*(5/9)
   
    #Converting pressure to bar
    pres=pres*0.0689476
    
    #Converting volumetric flowrate of liquid to m3/s
    vol_l=(bopd+bwpd)*1.84E-6
        
    #Converting volumetric flowrate of gas to m3/s 
    vol_g=(mscf*28.3168)/(24*3600)
    
    #Calculation of liquid and gas superficial velocity (m/s)
    v_sl=vol_l/((3.1416*diameter**2)/4)
    
    v_sg=(vol_g/((3.1416*diameter**2)/4))*0.9*((temp+273.15)/(288.7))
      
    #Calculation of superficial velocity of mixture
    v_m = v_sg + v_sl
    
    #Calculation of liquid fraction (holdup)
    holdup=vol_l/(vol_l+vol_g)
    
    #density of water [kg/m3]. Default value in Norsok: 1024 kg/m3
    den_w=998.2/(1+(0.0002*(temp-20)))
    
    #density of oil [kg/m3]. Default value in Norsok: 850 kg/m3
    den_o=850-((1.825-(0.001315*850))*(temp-20))
    
    #Calculation density of liquid
    water_cut=bwpd/(bwpd+bopd)
    density_l=(water_cut*den_w)+(den_o*(1-water_cut))
    
    #Calculation of density of gas. Default in Norsok: 10 kg/m3
    #density_g=28.97*((pres*0.75)/(0.9*8.314*(temp+273.25)*1000))
    density_g=(2.7*14.5*16.018*pres*0.5537)/(0.9*(460+(1.8*temp+32)))
    
    #Calculation of density of mixture
    density_mix = density_l * holdup + density_g * (1 - holdup)
    
    #Viscosity of gas conversion from cp to kg/m.s. I use the default value of Norsok (given in cp) for viscosity of gas. 
    vis_g=0.03/1000
    
    #Viscosity of oil from https://www.petroskills.com/blog/entry/crude-oil-and-changing-temperature#.X8i1prlxc2w for a crude of 25 API
    #Viscosity equation from "Estimating the Viscosity of Crude Oil System", SPE5434
    sg=0.9042-(5.93E-4)*(temp-15)
    api=(141.5/sg)-131.5
    y=10**(3.0324-0.02023*api)
    x=y*(1.8*temp+32)**-1.163
    vis_o=((10**x)-1)/1000
    
    #Viscosity of water
    vis_w=((temp+273.15)-225.4)**-1.637
    
    #Calculation viscosity of liquid    
    if water_cut < 0.5:
        bo1=water_cut/(1.187*(1-(1/7.06)**0.4))
        vis_l=vis_o*((1+((water_cut/bo1)/(1.187-(water_cut/bo1))))**2.5)
    else:
        bo2=(1-water_cut)/(1.187*(1-((vis_w/vis_o)/7.06)**0.4))
        vis_l=vis_w*(1+(((1-water_cut)/bo2)/(1.187-((1-water_cut)/bo2))))**2.5
    
    #Calculation viscosity of mixture   
    vis_m = vis_l * holdup + vis_g * (1- holdup)
            
    #Calculation friction factor. I use the default value of roughness [m] of Norsok
    roughness=50E-6
    friction = 0.001375 * (1 + ((20000 * roughness / diameter) + (10 ** 6 * (vis_m / (v_m * diameter * density_mix)))) ** 0.33)
    
    #Calculation of shear stress
    tempo = 0.5 * density_mix * friction * v_m ** 2
    return tempo


# Funcion para calcular la constante Kt utilizada en el modelo NORSOK

# In[13]:


def Kt(temp):   
    # Converting temp to C
    temp = (temp-32)*(5/9)
    
    temp_table = (5, 15, 20, 40, 60, 80, 90, 120, 150)
    value_table = (0.42, 1.59, 4.762, 8.927, 10.695, 9.949, 6.250, 7.770, 5.203)
    
    for i, temp_i in enumerate(temp_table):
        
        if temp <= temp_i:
            temp_lower = temp_table[i-1]
            temp_upper = temp_table[i]
            value_lower = value_table[i-1]
            value_upper = value_table[i]
            break

    return round(value_upper - (temp_upper - temp) * (value_upper - value_lower)/(temp_upper - temp_lower), 3)


# Funcion para calcular la velocidad de corrosion

# In[14]:


def calcNorsok(temp,pres,bopd,bwpd,mscf,co2fraction,HCO3,Cl,Na,K,Mg,Ca,Sr,Ba,SO4,cac,diameter,IC_efficiency,situacion):
    phcal = ph(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
    fphcal = fpH_Cal((temp-32)*(5/9), phcal)
    fy=FugacityofCO2(co2fraction, pres, temp, bopd, bwpd, mscf)
    kt = Kt(temp)
    shs = Shearstress(mscf, bopd, bwpd, diameter, temp, pres)

    #Calculo de la velocidad de corrosion segun Norsok, en mm/year
    nk = kt * fy ** 0.62 * (shs /19) ** (0.146 + 0.0324 * math.log10(fy)) * fphcal

    #Transformar a mpy y tomar en cuenta la eficiencia del inhibidor de corrosion
    corr_ic=nk*39.4*(100-IC_efficiency)/100

    #Asignar niveles de riesgo segun norma NACE
    if corr_ic < 1:
        corr_risk ='Bajo'
    if corr_ic>=1 and corr_ic<5:
        corr_risk ='Moderado'
    if corr_ic>=5 and corr_ic<10:
        corr_risk = "Alto"
    if corr_ic >= 10:
        corr_risk = 'Muy alto'

    return nk, corr_ic,corr_risk


# Funcion para calcular el indice de saturacion

# In[15]:


def calCalcite(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf, situacion): 
    
    calcite_si=indiceSaturacion(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
    
    solid=ptb(pres, temp, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
       
    #Asignar niveles de riesgo     
    if calcite_si <= 0.5:
            scale_risk ='Bajo'
    if calcite_si>0.5 and calcite_si<=1.5:
            scale_risk ='Moderado'
    if calcite_si>1.5 and calcite_si<=2.5:
            scale_risk ='Alto'
    if calcite_si > 2.5:
            scale_risk ='Muy alto'

    return calcite_si, solid, scale_risk


# Funcion para graficar el perfil de la velocidad de corrosion

# In[16]:


def grahpNorskok(temp,temp1,pres,pres1,bopd,bwpd,mscf,co2fraction,HCO3,Cl,Na,K,Mg,Ca,Sr,Ba,SO4,cac,diameter,IC_efficiency, h1,well_name):
  
  temp_array = np.linspace(temp, temp1, 100)
  press_array = np.linspace(pres, pres1, 100)
  depth_array = np.linspace(0,h1,100)
  aux = 0
  kt_df = []
  yt_df = []
  fy_df = []
  shearstress_df = []
  ph_df = []
  fph_df = []
  nk_df = []
  corr_profile_risk=[]
  
  for i in temp_array:
    auxkt = Kt(i)
    kt_df.append(auxkt)
    auxyt = ytgas(co2fraction, bopd, bwpd, mscf, i, press_array[aux]) 
    yt_df.append(auxyt)
    auxfy = FugacityofCO2(co2fraction, press_array[aux], i, bopd, bwpd, mscf)
    fy_df.append(auxfy)
    auxss = Shearstress(mscf, bopd, bwpd, diameter, i, press_array[aux])
    shearstress_df.append(auxss)
    auxph = ph(press_array[aux],i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
    ph_df.append(auxph)
    auxfph = fpH_Cal((i-32)*(5/9), auxph)
    fph_df.append(auxfph)
    auxnk = auxkt * auxfy ** 0.62 * (auxss /19) ** (0.146 + 0.0324 * math.log10(auxfy)) * auxfph*39.4*(100-IC_efficiency)/100
    nk_df.append(auxnk)
    
    #Asignar niveles de riesgo segun norma NACE
    if auxnk < 1:
        auxcorr ='Bajo'
    if auxnk>=1 and auxnk<5:
        auxcorr ='Moderado'
    if auxnk>=5 and auxnk<10:
        auxcorr = "Alto"
    if auxnk >= 10:
        auxcorr = 'Muy alto'
    
    corr_profile_risk.append(auxcorr)
    
    aux = aux + 1

  return temp_array,press_array,depth_array,fy_df,ph_df,nk_df,corr_profile_risk


# Funcion para graficar el perfil del indice de saturacion

# In[17]:


def graphCalcite(temp,temp1,pres,pres1,h1,bopd,bwpd,mscf,co2fraction,HCO3,Cl,Na,K,Mg,Ca,Sr,Ba,SO4,cac,well_name):
  temp_array = np.linspace(temp, temp1, 100)
  press_array = np.linspace(pres, pres1, 100)
  depth_array=np.linspace(0,h1,100)
  aux = 0
  alkalinity = []
  fy = []
  yt = []
  ph1 = []
  calcite = []
  ptb1 = []
  scale_profile_risk=[]

  for i in temp_array:
      auxak = hco(press_array[aux], i, Na, K, Mg, Ca, Sr,
                    Ba, Cl, SO4, HCO3, co2fraction, cac, bopd, bwpd, mscf)
      alkalinity.append(auxak)
      auxfy = fugacityCO2(i, press_array[aux])
      fy.append(auxfy)
      auxyt = ytgas(co2fraction, bopd, bwpd, mscf, i, press_array[aux])
      yt.append(auxyt)
      auxph = ph(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
      ph1.append(auxph)
      auxcalcite = indiceSaturacion(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
      calcite.append(auxcalcite)
      auxptb = ptb(press_array[aux], i, Na, K, Mg, Ca, Sr, Ba, Cl, SO4, HCO3, co2fraction, cac,bopd,bwpd,mscf)
      ptb1.append(auxptb)
        
      #Asignar niveles de riesgo     
      if auxcalcite <= 0.5:
            auxscale ='Bajo'
      if auxcalcite>0.5 and auxcalcite<=1.5:
            auxscale ='Moderado'
      if auxcalcite>1.5 and auxcalcite<=2.5:
            auxscale ='Alto'
      if auxcalcite > 2.5:
            auxscale ='Muy alto'
    
      scale_profile_risk.append(auxscale)
        
      aux = aux + 1
   
  return temp_array, press_array, depth_array, fy, ph1, calcite, ptb1, scale_profile_risk


# PRUEBA DE LAS FUNCIONES
# Modelo validado con la hoja de calculo Excel de Norsok 2017
# 
# Condiciones del pozo 
# Variables pres y temp corresponden a presion y temperatura en la cabeza del pozo (wellhead). 
# Variables pres1 y temp1 corresponden a presion y temperatura en fondo del pozo (bottom).

# Inicializar data frames para guardar los resultados 

# In[18]:


df0=[]
df1=[]
df2=[]
df3=[]
df4=[]
df5=[]
df6=[]
df7=[]
df8=[]
df9=[]
df10=[]
df11=[]
df12=[]
df13=[]
df14=[]
df15=[]
df16=[]
df17=[]
df18=[]
df19=[]
df20=[]
df21=[]
df22=[]
df23=[]
df24=[]
df25=[]
df26=[]
df27=[]
df28=[]
df29=[]
df30=[]
df31=[]
df32=[]
df33=[]
df34=[]
df35=[]
df36=[]
df37=[]
df38=[]
df39=[]
df40=[]
df41=[]
df42=[]
df43=[]
df44=[]
df45=[]
df46=[]
df47=[]

# Calculo de los resultados

# In[19]:


def run():
    
    st.set_page_config(layout="wide")
    
    from PIL import Image
    image = Image.open('logo_Pungo_1.png')
    
    col1,col2=st.columns(2)
    
    with col1:
        st.image(image)
    
    with col2: 
        st.header('DIGITALIZACIÓN QUE GENERA IMPACTO: Desarrollamos soluciones digitales que te permiten optimizar tus recursos (tiempo, dinero) y mejorar el desempeno de los pozos petroleros')
    
    add_selectbox = st.sidebar.selectbox(
    "Que deseas predecir?",
    ("Un solo pozo", "Varios pozos"))
    
    st.sidebar.info('Esta es una aplicación web para predecir la velocidad de corrosión, índice de saturacion, y criticidad de los pozos de petróleo en función de\
         varios parametros que puedes ver en la pantalla. Por favor, ajusta el\
         valor de cada parametro. Después de eso, da clic en el botón Calcular en la parte inferior para\
         ver las predicciones del modelo.')
    
    st.sidebar.info('Esta herramienta tambien permite optimizar las dosis de los quimicos anticorrosivo y antiescala en funcion de los resultados\
            obtenidos en el modelo (riesgo de corrosion, incrustaciones y criticidad). Para ello, da clic en el boton Optimizar en la parte inferior\
            y obtendras las dosis recomendadas de los quimicos, junto con el ahorro generado')
    
    st.title("Cálculo de velocidad de corrosión, índice de saturacion, y criticidad de pozos petroleros")

    if add_selectbox == 'Un solo pozo':
        
        st.subheader("Predicciones para un solo pozo")
             
        BOPD = st.slider(label = 'Barriles de petroleo por dia, BPPD', min_value = 0,
                          max_value = 5000 ,
                          value = 2500,step=1)
              
        BWPD = st.slider(label = 'Barriles de agua por dia, BWPD', min_value = 0,
                          max_value = 5000 ,
                          value = 2500,step=1)
                          
        MSCF = st.slider(label = 'Caudal de gas, MSCFD', min_value = 0.0,
                          max_value = 1000.0 ,
                          value = 500.0,step=0.1)
   
        pressure_head = st.slider(label = 'Presion de cabeza, psi', min_value = 0.0,
                          max_value = 600.0,
                          value = 150.0,step=0.1)

        temperature_head = st.slider(label = 'Temperatura de cabeza, F', min_value = 0.0,
                          max_value = 300.0 ,
                          value = 150.0,step=0.1)
        
        pressure_bottom = st.slider(label = 'Presion de fondo, psi', min_value = 0.0,
                          max_value = 5000.0,
                          value = 2500.0,step=0.1)

        temperature_bottom = st.slider(label = 'Temperatura de fondo, F', min_value = 0.0,
                          max_value = 300.0,
                          value = 250.0,step=0.1)

        chlorides = st.slider(label = 'Cloruros, ppm', min_value = 0,
                          max_value = 100000 ,
                          value = 50000,step=1)

        co2_gas = st.slider(label = 'CO2 gas, fraccion', min_value = 0.000,
                          max_value = 1.000 ,
                          value = 0.200,step=0.001)
                          
        alkalinity = st.slider(label = 'Alcalinidad, ppm', min_value = 0,
                          max_value = 5000,
                          value = 500,step=1)
        
        sodium = st.slider(label = 'Sodio, ppm', min_value = 0,
                          max_value = 50000,
                          value = 1000,step=1)
        
        potassium = st.slider(label = 'Potasio, ppm', min_value = 0,
                          max_value = 5000,
                          value = 500,step=1)
        
        magnesium = st.slider(label = 'Magnesio, ppm', min_value = 0,
                          max_value = 5000,
                          value = 500,step=1)
        
        calcium = st.slider(label = 'Calcio, ppm', min_value = 0,
                          max_value = 50000,
                          value = 500,step=1)
        
        strontium = st.slider(label = 'Estroncio, ppm', min_value = 0,
                          max_value = 5000,
                          value = 100,step=1)
        
        barium = st.slider(label = 'Bario, ppm', min_value = 0,
                          max_value = 5000,
                          value = 100,step=1)
        
        sulphates = st.slider(label = 'Sulfatos, ppm', min_value = 0,
                          max_value = 5000,
                          value = 100,step=1)
        
        carb_acids = st.slider(label = 'Acidos carboxilicos, ppm', min_value = 0,
                          max_value = 5000,
                          value = 100,step=1)
        
        pipe_diameter = st.slider(label = 'Diametro tuberia, in', min_value = 3.0,
                          max_value = 5.0,
                          value = 3.5,step=0.1)
     
        well_depth = st.slider(label = 'Profundidad pozo, ft', min_value = 0,
                          max_value = 20000,
                          value = 10000,step=1)
        
        dosis_ic = st.slider(label = 'Dosis Anticorrosivo, gal/dia', min_value = 0,
                          max_value = 20,
                          value = 10,step=1)
         
        dosis_is = st.slider(label = 'Dosis Antiescala, gal/dia', min_value = 0,
                          max_value = 20,
                          value = 10,step=1)
        
        #Asumo una eficiencia del inhibidor de corrosion del 97%
        IC_eff = 97
        
        i='Pozo_1'
        
        output=""

        features = {'BOPD': BOPD, 'BWPD': BWPD,
            'Caudal_gas_MSCFD': MSCF, 'Presion_cabeza_psi': pressure_head,
            'Temperatura_cabeza_F': temperature_head, 'Presion_fondo_psi': pressure_bottom,
            'Temperatura_fondo_F': temperature_bottom,'Cloruros_ppm': chlorides,
            'CO2_gas': co2_gas, 'Alcalinidad_ppm': alkalinity,
            'Sodio_ppm': sodium, 'Potasio_ppm': potassium,
            'Magnesio_ppm': magnesium, 'Calcio_ppm': calcium,
            'Estroncio_ppm': strontium, 'Bario_ppm': barium,
            'Sulfatos_ppm': sulphates, 'Acidos_carboxilicos_ppm': carb_acids,
            'Diametro tuberia_ppm': pipe_diameter, 'Eficiencia IC_%': IC_eff,
            'Profundidad pozo_ft':well_depth}
        
        features_df  = pd.DataFrame([features])
        
        #st.dataframe(features_df)
        
        #Velocidad de corrosion en cabeza
        nk_temp, corr_ic_temp,corr_risk_temp = calcNorsok(temperature_head,
           pressure_head,BOPD,BWPD,MSCF,co2_gas,alkalinity,chlorides,sodium,
           potassium,magnesium,calcium,strontium,barium,sulphates,carb_acids,
           pipe_diameter,IC_eff,'CABEZA')

        #Velocidad de corrosion en fondo
        nk_temp1, corr_ic_temp1,corr_risk_temp1 = calcNorsok(temperature_bottom,
           pressure_bottom,BOPD,BWPD,MSCF,co2_gas,alkalinity,chlorides,sodium,
           potassium,magnesium,calcium,strontium,barium,sulphates,carb_acids,
           pipe_diameter,IC_eff,'FONDO')

        #Indice de saturacion en cabeza
        calcite_si_temp, solid_temp,scale_risk_temp=calCalcite(pressure_head,
           temperature_head,sodium,potassium,magnesium,
           calcium,strontium,barium,chlorides,sulphates,
           alkalinity,co2_gas,carb_acids,BOPD,
           BWPD,MSCF,"CABEZA")

        #Indice de saturacion en fondo
        calcite_si_temp1, solid_temp1,scale_risk_temp1=calCalcite(pressure_bottom,
           temperature_bottom,sodium,potassium,magnesium,
           calcium,strontium,barium,chlorides,sulphates,
           alkalinity,co2_gas,carb_acids,BOPD,
           BWPD,MSCF,"FONDO")
        
        #Perfil de la velocidad de corrosion
        temp_array,press_array,depth_array,fy_df,ph_df,nk_df,corr_profile_risk = grahpNorskok(temperature_head,
           temperature_bottom,pressure_head,pressure_bottom,BOPD,BWPD,MSCF,co2_gas,alkalinity,
           chlorides,sodium,potassium,magnesium,calcium,strontium,barium,sulphates,carb_acids,
           pipe_diameter,IC_eff,well_depth,i)

        #Perfil del indice de saturacion
        temp_array, press_array, depth_array, fy, ph1, calcite, ptb1,scale_profile_risk=graphCalcite(temperature_head,
            temperature_bottom,pressure_head,pressure_bottom,
            well_depth,BOPD,BWPD,MSCF,co2_gas,alkalinity,
            chlorides,sodium,potassium,magnesium,calcium,
            strontium,barium,sulphates,carb_acids,i)
        
        #Guardar los resultados de cabeza y fondo en un data frame
        df0.append(i)
        df1.append(corr_ic_temp)
        df20.append(corr_risk_temp)
        df2.append(corr_ic_temp1)
        df21.append(corr_risk_temp1)
        df3.append(calcite_si_temp)
        df22.append(scale_risk_temp)
        df4.append(calcite_si_temp1)
        df23.append(scale_risk_temp1)
        
        #Resultados de velocidad de corrosion e indice de saturacion en cabeza y fondo
        results_corr=pd.DataFrame({'Velocidad de corrosion cabeza [mpy]':df1,
                                          'Riesgo de corrosion cabeza':df20,
                                          'Velocidad de corrosion fondo [mpy]':df2,
                                          'Riesgo de corrosion fondo':df21})
        
        results_esc=pd.DataFrame({'Indice de saturacion cabeza':df3,
                                          'Riesgo de incrustaciones cabeza':df22,
                                          'Indice de saturacion fondo':df4,
                                          'Riesgo de incrustaciones fondo':df23})  
        
        #Guardar los resultados del perfil de velocidad de corrosion en un data frame
        df10.append(temp_array)
        df11.append(press_array)
        df12.append(depth_array)
        df13.append(fy_df)
        df14.append(ph_df)
        df15.append(nk_df)
        df25.append(corr_profile_risk)
        
        #Resultados del perfil de velocidad de corrosion
        results_corr_profile=pd.DataFrame({'Pozo':df0,'Temperatura [F]':df10,'Presion [psi]':df11,
                              'Profundidad [ft]':df12,'Fugacidad CO2':df13,
                              'pH':df14,'Velocidad de corrosion (mpy)':df15,
                              'Riesgo de corrosion':df25}).set_index(['Pozo']).apply(pd.Series.explode).reset_index()  

        #Guardar los resultados del perfil del indice de saturacion en un data frame
        df16.append(fy)
        df17.append(ph1)
        df18.append(calcite)
        df19.append(ptb1)
        df26.append(scale_profile_risk)
        
        #Resultados del perfil del indice de saturacion
        results_scale_profile=pd.DataFrame({'Pozo':df0,'Temperatura [F]':df10,'Presion [psi]':df11,
                              'Profundidad [ft]':df12,'Fugacidad CO2':df16,
                              'pH':df17,'Indice de saturacion calcita':df18,'Solidos [PTB]':df19,
                              'Riesgo de incrustaciones':df26}).set_index(['Pozo']).apply(pd.Series.explode).reset_index()
        
        #Graficar el perfil de velocidad de corrosion
        fig_corr=px.line(results_corr_profile,x='Velocidad de corrosion (mpy)',y='Profundidad [ft]',title='Perfil de velocidad de corrosion',
            hover_name='Pozo',hover_data=['Presion [psi]','Temperatura [F]','Riesgo de corrosion'])

        fig_corr.update_traces(mode="markers+lines")
        fig_corr.update_xaxes(showspikes=True, spikecolor='black')
        fig_corr.update_yaxes(showspikes=True,spikecolor='black')
        fig_corr.update_yaxes(autorange="reversed")
        
        #Graficar el perfil del indice de saturacion
        fig_sca=px.line(results_scale_profile,x='Indice de saturacion calcita', y='Profundidad [ft]',title='Perfil del indice de saturacion',
                    hover_name='Pozo',hover_data=['Presion [psi]','Temperatura [F]','Solidos [PTB]','Riesgo de incrustaciones'])

        fig_sca.update_traces(mode="markers+lines")
        fig_sca.update_xaxes(showspikes=True, spikecolor='black')
        fig_sca.update_yaxes(showspikes=True,spikecolor='black')
        fig_sca.update_yaxes(autorange="reversed")
        
        #Calculo de la criticidad  
        prod=BOPD
        corr_median=np.median(df15)
        scale_median=np.median(df18)

        if prod > 500:
            critic_prod=4
        if prod >350 and prod<=500:
            critic_prod=3
        if prod >200 and prod<=350:
            critic_prod=2
        if prod <= 200:
            critic_prod=1

        if corr_median > 10:
            critic_corr=4
            risk_corr_median = 'Muy alto'
        if corr_median >5 and corr_median<=10:
            critic_corr=3
            risk_corr_median = "Alto"
        if corr_median >1 and corr_median<=5:
            critic_corr=2
            risk_corr_median ='Moderado'
        if corr_median <= 1:
            critic_corr=1
            risk_corr_median ='Bajo'

        if scale_median > 2.5:
            critic_si=4
            risk_scale_median = 'Muy alto'
        if scale_median >1.5 and scale_median<=2.5:
            critic_si=3
            risk_scale_median = 'Alto'
        if scale_median >0.5 and scale_median<=1.5:
            critic_si=2
            risk_scale_median = 'Moderado'
        if scale_median <= 0.5:
            critic_si=1
            risk_scale_median = 'Bajo'

        critic_tot=critic_prod*critic_corr*critic_si

        if critic_tot>32:
            nivel_critic='Muy alta'
        if critic_tot >16 and critic_tot<=32:
            nivel_critic='Alta'
        if critic_tot >8 and critic_tot<=16:
            nivel_critic='Moderada'
        if critic_tot<=8:
            nivel_critic='Baja'   

        #Guardar los resultados de la criticidad
        df5.append(BOPD)
        df6.append(critic_prod)
        df7.append(critic_corr)
        df8.append(critic_si)
        df9.append(critic_tot)
        df24.append(nivel_critic)
          
        
        #Calculo de dosis de quimicos y ahorro
        precio_ic=20
        precio_is=20

        if corr_median>10:
            dosis_recomendada_ic= math.ceil(80*BWPD/23810)
        if corr_median >5 and corr_median<=10:
            dosis_recomendada_ic= math.ceil(60*BWPD/23810)
        if corr_median >1 and corr_median<=5:
            dosis_recomendada_ic= math.ceil(40*BWPD/23810)
        if corr_median<=1:
            dosis_recomendada_ic= math.ceil(20*BWPD/23810) 

        diferencia_dosis_ic=dosis_ic - dosis_recomendada_ic

        if dosis_ic < dosis_recomendada_ic:
            estado_dosis_ic = "Sub-dosificado"
        else:
            estado_dosis_ic = "Sobre-dosificado"

        ahorro_anual_ic=diferencia_dosis_ic * precio_ic *30*12
        
        if ahorro_anual_ic>0:
            resultado_ic=str('Ahorro potencial '+"%.2f" % ahorro_anual_ic) + ' USD/año' + ' al disminuir dosis de anticorrosivo'
        else:
            resultado_ic='Pozo con riesgo de corrosion. Aumentar dosis de anticorrosivo'

        #Guardar los resultados de la dosis de anticorrosivo en data frames
        df27.append(dosis_ic)
        df28.append(dosis_recomendada_ic)
        df29.append(ahorro_anual_ic)
        df40.append(estado_dosis_ic)

        results_ic=pd.DataFrame({'Dosis actual de anticorrosivo [gal/dia]':df27,
                                 'Dosis recomendada de anticorrosivo [gal/dia]':df28,
                                 'Estado de la dosificacion de anticorrosivo':df40,
                                 'Ahorro por optimizacion de anticorrosivo [USD/año]':df29})


        if scale_median>2.5:
            dosis_recomendada_is= math.ceil(80*BWPD/23810)
        if scale_median >1.5 and scale_median<=2.5:
            dosis_recomendada_is= math.ceil(60*BWPD/23810)
        if scale_median >0.5 and scale_median<=1.5:
            dosis_recomendada_is= math.ceil(40*BWPD/23810)
        if scale_median<=0.5:
            dosis_recomendada_is= math.ceil(20*BWPD/23810) 

        diferencia_dosis_is=dosis_is - dosis_recomendada_is

        if dosis_is < dosis_recomendada_is:
            estado_dosis_is = "Sub-dosificado"
        else:
            estado_dosis_is = "Sobre-dosificado"

        ahorro_anual_is=diferencia_dosis_is * precio_is *30*12
        
        if ahorro_anual_is>0:
            resultado_is=str('Ahorro potencial '+"%.2f" % ahorro_anual_is) + ' USD/año' + ' al disminuir dosis de antiescala'
        else:
            resultado_is='Riesgo de incrustaciones. Aumentar dosis de antiescala'

        #Guardar los resultados de la dosis de antiescala en data frames
        df30.append(dosis_is)
        df31.append(dosis_recomendada_is)
        df32.append(ahorro_anual_is)
        df41.append(estado_dosis_is)

        results_is=pd.DataFrame({'Dosis actual de antiescala [gal/dia]':df30,
                                 'Dosis recomendada de antiescala [gal/dia]':df31,
                                 'Estado de la dosificacion de antiescala':df41,
                                 'Ahorro por optimizacion de antiescala [USD/año]':df32})
  
        ahorro_total=ahorro_anual_ic+ahorro_anual_is    

        if st.button('Calcular'):

            col1, col2, col3 = st.columns(3)
        
            with col1:
                st.metric('Velocidad de corrosion en cabeza',str("%.2f" % np.float_(df1[0]))+' MPY')
                st.metric('Riesgo de corrosion en cabeza',df20[0])

            with col2:
                st.metric('Velocidad de corrosion en fondo',str("%.2f" % np.float_(df2[0]))+' MPY')
                st.metric('Riesgo de corrosion en fondo',df21[0])

            with col3:
                st.metric('Velocidad de corrosion promedio',str("%.2f" % np.float_(corr_median))+' MPY')
                st.metric('Riesgo de corrosion promedio',risk_corr_median)

            st.plotly_chart(fig_corr, use_container_width=True)

            st.write("")    

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric('Indice de saturacion en cabeza',str("%.2f" % np.float_(df3[0])))
                st.metric('Riesgo de incrustaciones en cabeza',df22[0])

            with col2:
                st.metric('Indice de saturacion en fondo',str("%.2f" % np.float_(df4[0])))
                st.metric('Riesgo de incrustaciones en fondo',df23[0])

            with col3:
                st.metric('Indice de saturacion promedio',str("%.2f" % np.float_(scale_median)))
                st.metric('Riesgo de incrustaciones promedio',risk_scale_median)

            st.plotly_chart(fig_sca, use_container_width=True)

            col1,col2,col3=st.columns(3)

            with col1:
                st.metric('Produccion de petroleo del pozo',str("%.0f" % np.float_(BOPD))+' BOPD')

            with col2:
                st.metric('Indice de criticidad total del pozo',str("%.2f" % np.float_(df9[0])))

            with col3:
                st.metric('Prioridad para tratamiento quimico',nivel_critic)

            st.success('El indice de criticidad se calcula a partir de la produccion de petroleo del pozo, la velocidad de corrosion promedio, y el indice de\
                       saturacion promedio. Por tanto, es un valor que toma en cuenta el riesgo global de la perdida de produccion en caso de sufrir eventos de\
                       corrosion o incrustaciones a lo largo de toda la tuberia de produccion. El valor del indice de criticidad permite priorizar la dosificacion\
                       de quimicos y el monitoreo a los pozos mas criticos del campo. De esta manera, se optimiza el tiempo y los recursos del personal de\
                       tratamiento quimico.')
                
    
        if st.button('Optimizar dosis de quimicos'):

            #Dosis anticorrosivo
            col1,col2,col3=st.columns(3)

            with col1:
                st.metric('Dosis actual anticorrosivo',str("%.0f" % np.float_(dosis_ic))+' gal/dia')

            with col2:
                st.metric('Dosis recomendada anticorrosivo',str(dosis_recomendada_ic)+' gal/dia')

            with col3:
                st.metric('Estado inyeccion anticorrosivo',estado_dosis_ic)

            st.success(resultado_ic)

            #Dosis antiescala
            col1,col2,col3=st.columns(3)

            with col1:
                st.metric('Dosis actual antiescala',str("%.0f" % np.float_(dosis_is))+' gal/dia')

            with col2:
                st.metric('Dosis recomendada antiescala',str(dosis_recomendada_is)+' gal/dia')

            with col3:
                st.metric('Estado inyeccion antiescala',estado_dosis_is)

            st.success(resultado_is)

            st.markdown('Nota: se asume un precio de 20 USD/gal para el anticorrosivo y el antiescala')
     
    
        if st.button('Beneficios'): 
            st.success("Felicitaciones, has comenzado a desbloquear algunos de los beneficios de la digitalización y el procesamiento de datos con modelos cientificos e inteligencia artificial.") 
            st.markdown("En esta ocasion, te hemos presentado los beneficios del modulo de Calculos. Sin embargo, quedan inmensas oportunidades para optimizar la productividad y el desempeño de los pozos petroleros al integrar todos los modulos de ASTRO") 
            st.markdown("Si estás listo para emprender este camino y llevar la digitalización al siguiente nivel, contactate con nosotros para guiarte en el proceso y recibir asesoria de nuestros expertos para optimizar tus operaciones.")
    
            contact_form = """
            <form action="https://formsubmit.co/zapaz.consultores@gmail.com" method="POST">
                 <input type="hidden" name="_captcha" value="false">
                 <input type="text" name="Nombre" placeholder="Tu nombre" required>
                 <input type="email" name="email" placeholder="Tu email" required>
                 <textarea name="message" placeholder="Tu mensaje"></textarea>
                 <button type="submit">Send</button>
            </form>
            """

            st.markdown(contact_form, unsafe_allow_html=True)

            # Use Local CSS File
            def local_css(file_name):
                with open(file_name) as f:
                    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


            local_css("style.css")

    
    # Codigo para varios pozos
    
    if add_selectbox == 'Varios pozos':
        
        st.subheader("Predicciones para múltiples pozos")
        
        file_upload = st.file_uploader("Subir archivo CSV para las predicciones", type=["csv"])

        if file_upload is not None:
            data = pd.read_csv(file_upload)
            
            parameters = dict()
            parameters_names = np.array(["BOPD", "BWPD", "MSCF", "Temperature Head", "Temperature Bottom", "Pressure Head", "Pressure Bottom",
                    "CO2 Fraction in Gas", "Alkalinity", "Chlorides", "Sodium", "Magnesium", "Potassium","Calcium", "Strontium",
                    "Barium", "Sulphates", "Carboxylic acids", "Well Depth", "Well Pipe Diameter", "Corrosion Inhibitor Efficiency",
                    "aux","dosis_ic","dosis_is","precio_ic","precio_is"])
            
            for i in data.drop(["Well", "Unit"], axis = 1).columns:
                for id_j, j in enumerate(parameters_names):
                    parameters[j] = float(data[i].iloc[int(np.where(data["Well"] == j)[0])])

                #Velocidad de corrosion en cabeza
                nk_temp, corr_ic_temp,corr_risk_temp = calcNorsok(parameters["Temperature Head"],
                           parameters["Pressure Head"],
                           parameters["BOPD"],
                           parameters["BWPD"],
                           parameters["MSCF"],
                           parameters["CO2 Fraction in Gas"],
                           parameters["Alkalinity"],
                           parameters["Chlorides"],
                           parameters["Sodium"],
                           parameters["Potassium"],
                           parameters["Magnesium"],
                           parameters["Calcium"],
                           parameters["Strontium"],
                           parameters["Barium"],
                           parameters["Sulphates"],
                           parameters["Carboxylic acids"],
                           parameters["Well Pipe Diameter"],
                           parameters["Corrosion Inhibitor Efficiency"],
                           'CABEZA')

                #Velocidad de corrosion en fondo
                nk_temp1, corr_ic_temp1,corr_risk_temp1 = calcNorsok(parameters["Temperature Bottom"],
                           parameters["Pressure Bottom"],
                           parameters["BOPD"],
                           parameters["BWPD"],
                           parameters["MSCF"],
                           parameters["CO2 Fraction in Gas"],
                           parameters["Alkalinity"],
                           parameters["Chlorides"],
                           parameters["Sodium"],
                           parameters["Potassium"],
                           parameters["Magnesium"],
                           parameters["Calcium"],
                           parameters["Strontium"],
                           parameters["Barium"],
                           parameters["Sulphates"],
                           parameters["Carboxylic acids"],
                           parameters["Well Pipe Diameter"],
                           parameters["Corrosion Inhibitor Efficiency"],
                           'FONDO')

                #Perfil de la velocidad de corrosion
                temp_array,press_array,depth_array,fy_df,ph_df,nk_df,corr_profile_risk = grahpNorskok(parameters["Temperature Head"],
                           parameters["Temperature Bottom"],
                           parameters["Pressure Head"],
                           parameters["Pressure Bottom"],
                           parameters["BOPD"],
                           parameters["BWPD"],
                           parameters["MSCF"],
                           parameters["CO2 Fraction in Gas"],
                           parameters["Alkalinity"],
                           parameters["Chlorides"],
                           parameters["Sodium"],
                           parameters["Potassium"],
                           parameters["Magnesium"],
                           parameters["Calcium"],
                           parameters["Strontium"],
                           parameters["Barium"],
                           parameters["Sulphates"],
                           parameters["Carboxylic acids"],
                           parameters["Well Pipe Diameter"],
                           parameters["Corrosion Inhibitor Efficiency"],
                           parameters["Well Depth"],
                           i)

                #Indice de saturacion en cabeza
                calcite_si_temp, solid_temp,scale_risk_temp=calCalcite(parameters["Pressure Head"],
                           parameters["Temperature Head"],
                           parameters["Sodium"],
                           parameters["Potassium"],
                           parameters["Magnesium"],
                           parameters["Calcium"],
                           parameters["Strontium"],
                           parameters["Barium"],
                           parameters["Chlorides"],
                           parameters["Sulphates"],
                           parameters["Alkalinity"],
                           parameters["CO2 Fraction in Gas"],
                           parameters["Carboxylic acids"],
                           parameters["BOPD"],
                           parameters["BWPD"],
                           parameters["MSCF"],
                           "CABEZA")

                #Indice de saturacion en fondo
                calcite_si_temp1, solid_temp1,scale_risk_temp1=calCalcite(parameters["Pressure Bottom"],
                           parameters["Temperature Bottom"],
                           parameters["Sodium"],
                           parameters["Potassium"],
                           parameters["Magnesium"],
                           parameters["Calcium"],
                           parameters["Strontium"],
                           parameters["Barium"],
                           parameters["Chlorides"],
                           parameters["Sulphates"],
                           parameters["Alkalinity"],
                           parameters["CO2 Fraction in Gas"],
                           parameters["Carboxylic acids"],
                           parameters["BOPD"],
                           parameters["BWPD"],
                           parameters["MSCF"],
                           "FONDO")

                #Perfil del indice de saturacion
                temp_array, press_array, depth_array, fy, ph1, calcite, ptb1,scale_profile_risk=graphCalcite(parameters["Temperature Head"],
                            parameters["Temperature Bottom"],
                            parameters["Pressure Head"],
                            parameters["Pressure Bottom"],
                            parameters["Well Depth"],
                            parameters["BOPD"],
                            parameters["BWPD"],
                            parameters["MSCF"],
                            parameters["CO2 Fraction in Gas"],
                            parameters["Alkalinity"],
                            parameters["Chlorides"],
                            parameters["Sodium"],
                            parameters["Potassium"],
                            parameters["Magnesium"],
                            parameters["Calcium"],
                            parameters["Strontium"],
                            parameters["Barium"],
                            parameters["Sulphates"],
                            parameters["Carboxylic acids"],
                            i)    
                
                #Guardar los resultados de cabeza y fondo en un data frame
                df0.append(i)
                df1.append(corr_ic_temp)
                df20.append(corr_risk_temp)
                df2.append(corr_ic_temp1)
                df21.append(corr_risk_temp1)
                df3.append(calcite_si_temp)
                df22.append(scale_risk_temp)
                df4.append(calcite_si_temp1)
                df23.append(scale_risk_temp1)
                df5.append(parameters['BOPD'])
                                     
                
                #Guardar los resultados del perfil de velocidad de corrosion en un data frame
                df10.append(temp_array)
                df11.append(press_array)
                df12.append(depth_array)
                df13.append(fy_df)
                df14.append(ph_df)
                df15.append(nk_df)
                df25.append(corr_profile_risk)
                
                #Guardar los resultados del perfil del indice de saturacion en un data frame
                df16.append(fy)
                df17.append(ph1)
                df18.append(calcite)
                df19.append(ptb1)
                df26.append(scale_profile_risk)                
                
                #Calculo de la criticidad
                
                prod=parameters['BOPD']
                
                #Velocidad de corrosion promedio
                for list in df15:
                    corr_median=np.median(list)
               
               #Indice de saturacion promedio
                for list in df18:
                    scale_median=np.median(list)
                

                if prod > 500:
                    critic_prod=4
                if prod >350 and prod<=500:
                    critic_prod=3
                if prod >200 and prod<=350:
                    critic_prod=2
                if prod <= 200:
                    critic_prod=1

                if corr_median > 10:
                    critic_corr=4
                    risk_corr_median = 'Muy alto'
                if corr_median >5 and corr_median<=10:
                    critic_corr=3
                    risk_corr_median = "Alto"
                if corr_median >1 and corr_median<=5:
                    critic_corr=2
                    risk_corr_median ='Moderado'
                if corr_median <= 1:
                    critic_corr=1
                    risk_corr_median ='Bajo'

                if scale_median > 2.5:
                    critic_si=4
                    risk_scale_median = 'Muy alto'
                if scale_median >1.5 and scale_median<=2.5:
                    critic_si=3
                    risk_scale_median = 'Alto'
                if scale_median >0.5 and scale_median<=1.5:
                    critic_si=2
                    risk_scale_median = 'Moderado'
                if scale_median <= 0.5:
                    critic_si=1
                    risk_scale_median = 'Bajo'

                critic_tot=critic_prod*critic_corr*critic_si

                if critic_tot>32:
                    nivel_critic='Muy alta'
                if critic_tot >16 and critic_tot<=32:
                    nivel_critic='Alta'
                if critic_tot >8 and critic_tot<=16:
                    nivel_critic='Moderada'
                if critic_tot<=8:
                    nivel_critic='Baja'
                
                #Guardar los resultados de la criticidad
                df6.append(critic_prod)
                df7.append(critic_corr)
                df8.append(critic_si)
                df9.append(critic_tot)
                df24.append(nivel_critic)
                df44.append(corr_median)
                df45.append(scale_median)
                df46.append(risk_corr_median)
                df47.append(risk_scale_median)
                
                #Calculo de la dosis de quimico recomendada y el ahorro      
                BWPD=parameters['BWPD']
                dosis_ic=parameters["dosis_ic"]
                dosis_is=parameters["dosis_is"]
                precio_ic=parameters["precio_ic"]
                precio_is=parameters["precio_is"]
                
                if corr_median>10:
                    dosis_recomendada_ic= math.ceil(80*BWPD/23810)
                if corr_median >5 and corr_median<=10:
                    dosis_recomendada_ic= math.ceil(60*BWPD/23810)
                if corr_median >1 and corr_median<=5:
                    dosis_recomendada_ic= math.ceil(40*BWPD/23810)
                if corr_median<=1:
                    dosis_recomendada_ic= math.ceil(20*BWPD/23810) 

                diferencia_dosis_ic=dosis_ic - dosis_recomendada_ic

                if dosis_ic < dosis_recomendada_ic:
                    estado_dosis_ic = "Sub-dosificado"
                else:
                    estado_dosis_ic = "Sobre-dosificado"

                ahorro_anual_ic=diferencia_dosis_ic * precio_ic *30*12

                if ahorro_anual_ic>0:
                    resultado_ic=str('Ahorro potencial '+"%.2f" % ahorro_anual_ic) + ' USD/año' + ' al disminuir dosis de anticorrosivo'
                else:
                    resultado_ic='Pozo con riesgo de corrosion. Aumentar dosis de anticorrosivo'
                
                
                if scale_median>2.5:
                    dosis_recomendada_is= math.ceil(80*BWPD/23810)
                if scale_median >1.5 and scale_median<=2.5:
                    dosis_recomendada_is= math.ceil(60*BWPD/23810)
                if scale_median >0.5 and scale_median<=1.5:
                    dosis_recomendada_is= math.ceil(40*BWPD/23810)
                if scale_median<=0.5:
                    dosis_recomendada_is= math.ceil(20*BWPD/23810) 

                diferencia_dosis_is=dosis_is - dosis_recomendada_is

                if dosis_is < dosis_recomendada_is:
                    estado_dosis_is = "Sub-dosificado"
                else:
                    estado_dosis_is = "Sobre-dosificado"

                ahorro_anual_is=diferencia_dosis_is * precio_is *30*12

                if ahorro_anual_is>0:
                    resultado_is=str('Ahorro potencial '+"%.2f" % ahorro_anual_is) + ' USD/año' + ' al disminuir dosis de antiescala'
                else:
                    resultado_is='Riesgo de incrustaciones. Aumentar dosis de antiescala'
                
                #Guardar los resultados de dosis recomendadas y ahorro en un data frame
                df33.append(dosis_ic)
                df34.append(dosis_recomendada_ic)
                df35.append(resultado_ic)
                df36.append(dosis_is)
                df37.append(dosis_recomendada_is)
                df38.append(resultado_is)
                df42.append(estado_dosis_ic)
                df43.append(estado_dosis_is)
                
                
                #Resultados de corrosion                

                results_corr=pd.DataFrame({'Pozo':df0,'Velocidad de corrosion cabeza [mpy]':df1,'Riesgo de corrosion cabeza':df20,
                                      'Velocidad de corrosion fondo [mpy]':df2,'Riesgo de corrosion fondo':df21,
                                      'Velocidad de corrosion promedio [mpy]':df44,'Riesgo de corrosion promedio':df46}).set_index(['Pozo'])

                results_corr_profile=pd.DataFrame({'Pozo':df0,'Temperatura [F]':df10,'Presion [psi]':df11,
                                      'Profundidad [ft]':df12,'Fugacidad CO2':df13,
                                      'pH':df14,'Velocidad de corrosion (mpy)':df15,
                                      'Riesgo de corrosion':df25}).set_index(['Pozo']).apply(pd.Series.explode).reset_index()

                for i, (Pozo, subdf) in enumerate(results_corr_profile.groupby('Pozo'), 1):
                    locals()[f'well_corr{i}'] = subdf
      

                #Resultados de escala
                
                results_scale=pd.DataFrame({'Pozo':df0,'Indice de saturacion cabeza':df3,'Riesgo de incrustaciones cabeza':df22,
                                      'Indice de saturacion fondo':df4,'Riesgo de incrustaciones fondo':df23,
                                      'Indice de saturacion promedio [mpy]':df45,'Riesgo de incrustaciones promedio':df47}).set_index(['Pozo'])

                results_scale_profile=pd.DataFrame({'Pozo':df0,'Temperatura [F]':df10,'Presion [psi]':df11,
                                      'Profundidad [ft]':df12,'Fugacidad CO2':df16,
                                      'pH':df17,'Indice de saturacion calcita':df18,'Solidos [PTB]':df19,
                                      'Riesgo de incrustaciones':df26}).set_index(['Pozo']).apply(pd.Series.explode).reset_index()

                for i, (Pozo, subdf) in enumerate(results_scale_profile.groupby('Pozo'), 1):
                    locals()[f'well_scale{i}'] = subdf
                  
                
                
                #Resultados criticidad
                
                results_critic=pd.DataFrame({'Pozo':df0,'Producción [bopd]':df5,'Velocidad de corrosion promedio [mpy]':df44,
                                            'Indice de saturacion promedio [mpy]':df45, 
                                            'Criticidad total':df9, 'Prioridad TQ':df24})
                   
                
              
                #Resultados optimizacion quimicos
                
                results_opt_corr=pd.DataFrame({'Pozo':df0,'Dosis actual de anticorrosivo [gal/dia]':df33,
                                         'Dosis recomendada de anticorrosivo [gal/dia]':df34,
                                         'Estado dosificacion de anticorrosivo':df42,
                                         'Analisis inyeccion de anticorrosivo':df35})            
                                 
                results_opt_scale=pd.DataFrame({'Pozo':df0,'Dosis actual de antiescala [gal/dia]':df36,
                                         'Dosis recomendada de antiescala [gal/dia]':df37,
                                         'Estado dosificacion de antiescala':df43, 
                                         'Analisis inyeccion de antiescala':df38})
                
                    
            if st.button('Resultados de corrosion'):
                
                st.dataframe(results_corr)
            
                def convert_df(df):
                    return df.to_csv().encode('utf-8')
                
                csv_corr = convert_df(results_corr)
                st.download_button("📥Press to Download",csv_corr,"file.csv","text/csv",key='download-csv')
                
                corr_sliced=[v for k, v in results_corr_profile.groupby('Pozo')]

                for df in corr_sliced:
                    fig_corr=px.line(df,x='Velocidad de corrosion (mpy)',y='Profundidad [ft]',title='Perfil de velocidad de corrosion',
                    hover_name='Pozo',hover_data=['Presion [psi]','Temperatura [F]','Riesgo de corrosion'])

                    fig_corr.update_traces(mode="markers+lines")
                    fig_corr.update_xaxes(showspikes=True, spikecolor='black')
                    fig_corr.update_yaxes(showspikes=True,spikecolor='black')
                    fig_corr.update_yaxes(autorange="reversed")
                    st.plotly_chart(fig_corr, use_container_width=True)
                
            if st.button('Resultados indice de saturacion'):
                
                st.dataframe(results_scale)
            
                def convert_df(df):
                    return df.to_csv().encode('utf-8')
                
                csv_scale = convert_df(results_scale)                
                st.download_button("📥Press to Download",csv_scale,"file.csv","text/csv",key='download-csv')
                
                scale_sliced=[v for k, v in results_scale_profile.groupby('Pozo')]

                for df in scale_sliced:
                    fig_sca=px.line(df,x='Indice de saturacion calcita',y='Profundidad [ft]',title='Perfil del indice de saturacion',
                    hover_name='Pozo',hover_data=['Presion [psi]','Temperatura [F]','Solidos [PTB]','Riesgo de incrustaciones'])

                    fig_sca.update_traces(mode="markers+lines")
                    fig_sca.update_xaxes(showspikes=True, spikecolor='black')
                    fig_sca.update_yaxes(showspikes=True,spikecolor='black')
                    fig_sca.update_yaxes(autorange="reversed")
                    st.plotly_chart(fig_sca, use_container_width=True)   
            
            if st.button('Criticidad de pozos'): 
                
                st.dataframe(results_critic)
            
                def convert_df(df):
                    return df.to_csv().encode('utf-8')
                
                csv_critic = convert_df(results_critic)
                st.download_button("📥Press to Download",csv_critic,"file.csv","text/csv",key='download-csv')
                
                fig_crit=px.bar(results_critic, x='Pozo', y='Criticidad total',hover_data=['Prioridad TQ'])
                fig_crit.update_layout(xaxis={'categoryorder':'total descending'}, title='Criticidad total Pozos')
                st.plotly_chart(fig_crit, use_container_width=True)
                
                y1=np.float_(df6)
                y2=np.float_(df7)
                y3=np.float_(df8)
                ytot=y1+y2+y3

                fig_con=go.Figure()
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y1/ytot)*100,name='Produccion'))
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y2/ytot)*100,name='Corrosion'))
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y3/ytot)*100,name='Escala'))
                fig_con.update_layout(barmode='stack', title='Contribuciones (%) a la Criticidad total',yaxis_title='Contribucion, %')
                st.plotly_chart(fig_con, use_container_width=True)
                
                
            if st.button('Optimizar dosis de quimicos'):
                
                st.dataframe(results_opt_corr)
                def convert_df(df):
                    return df.to_csv().encode('utf-8')
                
                csv_opt_corr = convert_df(results_opt_corr)
                st.download_button("📥Press to Download",csv_opt_corr,"file.csv","text/csv",key='download-csv')
                
                st.dataframe(results_opt_scale)
                def convert_df(df):
                    return df.to_csv().encode('utf-8')
                
                csv_opt_scale = convert_df(results_opt_scale)
                st.download_button("📥Press to Download",csv_opt_scale,"file.csv","text/csv",key='download-csv')
                
                #ahorro_total_opt=results_opt['Ahorro total [USD/año]'].sum()
                #output4=str("%.2f" % ahorro_total_opt) + ' USD/año'
                #st.success('El ahorro total por optimizacion de quimicos es {}'.format(output4))

                #st.info('Nota: se asume un precio de 20 USD/gal para el anticorrosivo y el antiescala')
            

# In[20]:


if __name__ == '__main__':
    run()


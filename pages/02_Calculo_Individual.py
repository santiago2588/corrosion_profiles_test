import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from Functions.funciones import *

# PRUEBA DE LAS FUNCIONES
# Modelo validado con la hoja de calculo Excel de Norsok 2017
# Condiciones del pozo
# Variables pres y temp corresponden a presion y temperatura en la cabeza del pozo (wellhead). 
# Variables pres1 y temp1 corresponden a presion y temperatura en fondo del pozo (bottom).

# Inicializar data frames para guardar los resultados

# In[18]:

df0 = []
df1 = []
df2 = []
df3 = []
df4 = []
df5 = []
df6 = []
df7 = []
df8 = []
df9 = []
df10 = []
df11 = []
df12 = []
df13 = []
df14 = []
df15 = []
df16 = []
df17 = []
df18 = []
df19 = []
df20 = []
df21 = []
df22 = []
df23 = []
df24 = []
df25 = []
df26 = []
df27 = []
df28 = []
df29 = []
df30 = []
df31 = []
df32 = []
df33 = []
df34 = []
df35 = []
df36 = []
df37 = []
df38 = []
df39 = []
df40 = []
df41 = []
df42 = []
df43 = []
df44 = []
df45 = []
df46 = []
df47 = []


# Calculo de los resultados

# In[19]:
def run():
    st.set_page_config(layout="wide")

    st.title('Plan Gratuito')
    st.subheader("Cálculo de velocidad de corrosión, índice de saturacion, y criticidad para un pozo petrolero")

    st.write('Por favor, sigue los pasos que se presentan a continuacion.')
    st.write("1. Ajusta el valor de cada parametro en las unidades y rangos correspondientes")
    st.write("2. Presiona el boton Calcular para ver los resultados: velocidad de corrosion, indice de saturacion y criticidad del pozo")
    st.write("3. Presiona el boton Optimizar para obtener las dosis recomendadas de los quimicos, junto con el ahorro generado")

    BOPD = st.slider(label='Barriles de petroleo por dia, BPPD', min_value=0,
                     max_value=5000,
                     value=2500, step=1)

    BWPD = st.slider(label='Barriles de agua por dia, BWPD', min_value=0,
                     max_value=5000,
                     value=2500, step=1)

    MSCF = st.slider(label='Caudal de gas, MSCFD', min_value=0.0,
                     max_value=1000.0,
                     value=500.0, step=0.1)

    pressure_head = st.slider(label='Presion de cabeza, psi', min_value=0.0,
                              max_value=600.0,
                              value=150.0, step=0.1)

    temperature_head = st.slider(label='Temperatura de cabeza, F', min_value=0.0,
                                 max_value=300.0,
                                 value=150.0, step=0.1)

    pressure_bottom = st.slider(label='Presion de fondo, psi', min_value=0.0,
                                max_value=5000.0,
                                value=2500.0, step=0.1)

    temperature_bottom = st.slider(label='Temperatura de fondo, F', min_value=0.0,
                                   max_value=300.0,
                                   value=250.0, step=0.1)

    chlorides = st.slider(label='Cloruros, ppm', min_value=0,
                          max_value=100000,
                          value=50000, step=1)

    co2_gas = st.slider(label='CO2 gas, fraccion', min_value=0.000,
                        max_value=1.000,
                        value=0.200, step=0.001)

    alkalinity = st.slider(label='Alcalinidad, ppm', min_value=0,
                           max_value=5000,
                           value=500, step=1)

    sodium = st.slider(label='Sodio, ppm', min_value=0,
                       max_value=50000,
                       value=1000, step=1)

    potassium = st.slider(label='Potasio, ppm', min_value=0,
                          max_value=5000,
                          value=500, step=1)

    magnesium = st.slider(label='Magnesio, ppm', min_value=0,
                          max_value=5000,
                          value=500, step=1)

    calcium = st.slider(label='Calcio, ppm', min_value=0,
                        max_value=50000,
                        value=500, step=1)

    strontium = st.slider(label='Estroncio, ppm', min_value=0,
                          max_value=5000,
                          value=100, step=1)

    barium = st.slider(label='Bario, ppm', min_value=0,
                       max_value=5000,
                       value=100, step=1)

    sulphates = st.slider(label='Sulfatos, ppm', min_value=0,
                          max_value=5000,
                          value=100, step=1)

    carb_acids = st.slider(label='Acidos carboxilicos, ppm', min_value=0,
                           max_value=5000,
                           value=100, step=1)

    pipe_diameter = st.slider(label='Diametro tuberia, in', min_value=3.0,
                              max_value=5.0,
                              value=3.5, step=0.1)

    well_depth = st.slider(label='Profundidad pozo, ft', min_value=0,
                           max_value=20000,
                           value=10000, step=1)

    dosis_ic = st.slider(label='Dosis Anticorrosivo, gal/dia', min_value=0,
                         max_value=20,
                         value=10, step=1)

    dosis_is = st.slider(label='Dosis Antiescala, gal/dia', min_value=0,
                         max_value=20,
                         value=10, step=1)

    # Asumo una eficiencia del inhibidor de corrosion del 97%
    IC_eff = 97

    i = 'Pozo_1'

    output = ""

    features = {'BOPD': BOPD, 'BWPD': BWPD,
                'Caudal_gas_MSCFD': MSCF, 'Presion_cabeza_psi': pressure_head,
                'Temperatura_cabeza_F': temperature_head, 'Presion_fondo_psi': pressure_bottom,
                'Temperatura_fondo_F': temperature_bottom, 'Cloruros_ppm': chlorides,
                'CO2_gas': co2_gas, 'Alcalinidad_ppm': alkalinity,
                'Sodio_ppm': sodium, 'Potasio_ppm': potassium,
                'Magnesio_ppm': magnesium, 'Calcio_ppm': calcium,
                'Estroncio_ppm': strontium, 'Bario_ppm': barium,
                'Sulfatos_ppm': sulphates, 'Acidos_carboxilicos_ppm': carb_acids,
                'Diametro tuberia_ppm': pipe_diameter, 'Eficiencia IC_%': IC_eff,
                'Profundidad pozo_ft': well_depth}

    features_df = pd.DataFrame([features])

    # st.dataframe(features_df)

    # Velocidad de corrosion en cabeza
    nk_temp, corr_ic_temp, corr_risk_temp = calcNorsok(temperature_head,
                                                       pressure_head, BOPD, BWPD, MSCF, co2_gas, alkalinity,
                                                       chlorides, sodium,
                                                       potassium, magnesium, calcium, strontium, barium, sulphates,
                                                       carb_acids,
                                                       pipe_diameter, IC_eff, 'CABEZA')

    # Velocidad de corrosion en fondo
    nk_temp1, corr_ic_temp1, corr_risk_temp1 = calcNorsok(temperature_bottom,
                                                          pressure_bottom, BOPD, BWPD, MSCF, co2_gas, alkalinity,
                                                          chlorides, sodium,
                                                          potassium, magnesium, calcium, strontium, barium,
                                                          sulphates, carb_acids,
                                                          pipe_diameter, IC_eff, 'FONDO')

    # Indice de saturacion en cabeza
    calcite_si_temp, solid_temp, scale_risk_temp = calCalcite(pressure_head,
                                                              temperature_head, sodium, potassium, magnesium,
                                                              calcium, strontium, barium, chlorides, sulphates,
                                                              alkalinity, co2_gas, carb_acids, BOPD,
                                                              BWPD, MSCF, "CABEZA")

    # Indice de saturacion en fondo
    calcite_si_temp1, solid_temp1, scale_risk_temp1 = calCalcite(pressure_bottom,
                                                                 temperature_bottom, sodium, potassium, magnesium,
                                                                 calcium, strontium, barium, chlorides, sulphates,
                                                                 alkalinity, co2_gas, carb_acids, BOPD,
                                                                 BWPD, MSCF, "FONDO")

    # Perfil de la velocidad de corrosion
    temp_array, press_array, depth_array, fy_df, ph_df, nk_df, corr_profile_risk = grahpNorskok(temperature_head,
                                                                                                temperature_bottom,
                                                                                                pressure_head,
                                                                                                pressure_bottom,
                                                                                                BOPD, BWPD, MSCF,
                                                                                                co2_gas, alkalinity,
                                                                                                chlorides, sodium,
                                                                                                potassium,
                                                                                                magnesium, calcium,
                                                                                                strontium, barium,
                                                                                                sulphates,
                                                                                                carb_acids,
                                                                                                pipe_diameter,
                                                                                                IC_eff, well_depth,
                                                                                                i)

    # Perfil del indice de saturacion
    temp_array, press_array, depth_array, fy, ph1, calcite, ptb1, scale_profile_risk = graphCalcite(
        temperature_head,
        temperature_bottom, pressure_head, pressure_bottom,
        well_depth, BOPD, BWPD, MSCF, co2_gas, alkalinity,
        chlorides, sodium, potassium, magnesium, calcium,
        strontium, barium, sulphates, carb_acids, i)

    # Guardar los resultados de cabeza y fondo en un data frame
    df0.append(i)
    df1.append(corr_ic_temp)
    df20.append(corr_risk_temp)
    df2.append(corr_ic_temp1)
    df21.append(corr_risk_temp1)
    df3.append(calcite_si_temp)
    df22.append(scale_risk_temp)
    df4.append(calcite_si_temp1)
    df23.append(scale_risk_temp1)

    # Resultados de velocidad de corrosion e indice de saturacion en cabeza y fondo
    results_corr = pd.DataFrame({'Velocidad de corrosion cabeza [mpy]': df1,
                                 'Riesgo de corrosion cabeza': df20,
                                 'Velocidad de corrosion fondo [mpy]': df2,
                                 'Riesgo de corrosion fondo': df21})

    results_esc = pd.DataFrame({'Indice de saturacion cabeza': df3,
                                'Riesgo de incrustaciones cabeza': df22,
                                'Indice de saturacion fondo': df4,
                                'Riesgo de incrustaciones fondo': df23})

    # Guardar los resultados del perfil de velocidad de corrosion en un data frame
    df10.append(temp_array)
    df11.append(press_array)
    df12.append(depth_array)
    df13.append(fy_df)
    df14.append(ph_df)
    df15.append(nk_df)
    df25.append(corr_profile_risk)

    # Resultados del perfil de velocidad de corrosion
    results_corr_profile = pd.DataFrame({'Pozo': df0, 'Temperatura [F]': df10, 'Presion [psi]': df11,
                                         'Profundidad [ft]': df12, 'Fugacidad CO2': df13,
                                         'pH': df14, 'Velocidad de corrosion (mpy)': df15,
                                         'Riesgo de corrosion': df25}).set_index(['Pozo']).apply(
        pd.Series.explode).reset_index()

    # Guardar los resultados del perfil del indice de saturacion en un data frame
    df16.append(fy)
    df17.append(ph1)
    df18.append(calcite)
    df19.append(ptb1)
    df26.append(scale_profile_risk)

    # Resultados del perfil del indice de saturacion
    results_scale_profile = pd.DataFrame({'Pozo': df0, 'Temperatura [F]': df10, 'Presion [psi]': df11,
                                          'Profundidad [ft]': df12, 'Fugacidad CO2': df16,
                                          'pH': df17, 'Indice de saturacion calcita': df18, 'Solidos [PTB]': df19,
                                          'Riesgo de incrustaciones': df26}).set_index(['Pozo']).apply(
        pd.Series.explode).reset_index()

    # Graficar el perfil de velocidad de corrosion
    fig_corr = px.line(results_corr_profile, x='Velocidad de corrosion (mpy)', y='Profundidad [ft]',
                       title='Perfil de velocidad de corrosion',
                       hover_name='Pozo', hover_data=['Presion [psi]', 'Temperatura [F]', 'Riesgo de corrosion'])

    fig_corr.update_traces(mode="markers+lines")
    fig_corr.update_xaxes(showspikes=True, spikecolor='black')
    fig_corr.update_yaxes(showspikes=True, spikecolor='black')
    fig_corr.update_yaxes(autorange="reversed")

    # Graficar el perfil del indice de saturacion
    fig_sca = px.line(results_scale_profile, x='Indice de saturacion calcita', y='Profundidad [ft]',
                      title='Perfil del indice de saturacion',
                      hover_name='Pozo',
                      hover_data=['Presion [psi]', 'Temperatura [F]', 'Solidos [PTB]', 'Riesgo de incrustaciones'])

    fig_sca.update_traces(mode="markers+lines")
    fig_sca.update_xaxes(showspikes=True, spikecolor='black')
    fig_sca.update_yaxes(showspikes=True, spikecolor='black')
    fig_sca.update_yaxes(autorange="reversed")

    # Calculo de la criticidad
    prod = BOPD
    corr_max = np.max(df15)
    scale_max = np.max(df18)

    if prod > 500:
        critic_prod = 4
    if prod > 350 and prod <= 500:
        critic_prod = 3
    if prod > 200 and prod <= 350:
        critic_prod = 2
    if prod <= 200:
        critic_prod = 1

    if corr_max > 10:
        critic_corr = 4
        risk_corr_max = 'Muy alto'
    if corr_max > 5 and corr_max <= 10:
        critic_corr = 3
        risk_corr_max = "Alto"
    if corr_max > 1 and corr_max <= 5:
        critic_corr = 2
        risk_corr_max = 'Moderado'
    if corr_max <= 1:
        critic_corr = 1
        risk_corr_max = 'Bajo'

    if scale_max > 2.5:
        critic_si = 4
        risk_scale_max = 'Muy alto'
    if scale_max > 1.5 and scale_max <= 2.5:
        critic_si = 3
        risk_scale_max = 'Alto'
    if scale_max > 0.5 and scale_max <= 1.5:
        critic_si = 2
        risk_scale_max = 'Moderado'
    if scale_max <= 0.5:
        critic_si = 1
        risk_scale_max = 'Bajo'

    critic_tot = critic_prod * critic_corr * critic_si

    if critic_tot > 32:
        nivel_critic = 'Muy alta'
    if critic_tot > 16 and critic_tot <= 32:
        nivel_critic = 'Alta'
    if critic_tot > 8 and critic_tot <= 16:
        nivel_critic = 'Moderada'
    if critic_tot <= 8:
        nivel_critic = 'Baja'

        # Guardar los resultados de la criticidad
    df5.append(BOPD)
    df6.append(critic_prod)
    df7.append(critic_corr)
    df8.append(critic_si)
    df9.append(critic_tot)
    df24.append(nivel_critic)

    # Calculo de dosis de quimicos y ahorro
    precio_ic = 20
    precio_is = 20

    if corr_max > 10:
        dosis_recomendada_ic = math.ceil(80 * BWPD / 23810)
    if corr_max > 5 and corr_max <= 10:
        dosis_recomendada_ic = math.ceil(60 * BWPD / 23810)
    if corr_max > 1 and corr_max <= 5:
        dosis_recomendada_ic = math.ceil(40 * BWPD / 23810)
    if corr_max <= 1:
        dosis_recomendada_ic = math.ceil(20 * BWPD / 23810)

    diferencia_dosis_ic = dosis_ic - dosis_recomendada_ic

    if dosis_ic < dosis_recomendada_ic:
        estado_dosis_ic = "Sub-dosificado"
    else:
        estado_dosis_ic = "Sobre-dosificado"

    ahorro_anual_ic = diferencia_dosis_ic * precio_ic * 30 * 12

    if ahorro_anual_ic > 0:
        resultado_ic = str(
            'Ahorro potencial ' + "%.2f" % ahorro_anual_ic) + ' USD/año' + ' al disminuir dosis de anticorrosivo'
    else:
        resultado_ic = 'Pozo con riesgo de corrosion. Aumentar dosis de anticorrosivo'

    # Guardar los resultados de la dosis de anticorrosivo en data frames
    df27.append(dosis_ic)
    df28.append(dosis_recomendada_ic)
    df29.append(ahorro_anual_ic)
    df40.append(estado_dosis_ic)

    results_ic = pd.DataFrame({'Dosis actual de anticorrosivo [gal/dia]': df27,
                               'Dosis recomendada de anticorrosivo [gal/dia]': df28,
                               'Estado de la dosificacion de anticorrosivo': df40,
                               'Ahorro por optimizacion de anticorrosivo [USD/año]': df29})

    if scale_max > 2.5:
        dosis_recomendada_is = math.ceil(80 * BWPD / 23810)
    if scale_max > 1.5 and scale_max <= 2.5:
        dosis_recomendada_is = math.ceil(60 * BWPD / 23810)
    if scale_max > 0.5 and scale_max <= 1.5:
        dosis_recomendada_is = math.ceil(40 * BWPD / 23810)
    if scale_max <= 0.5:
        dosis_recomendada_is = math.ceil(20 * BWPD / 23810)

    diferencia_dosis_is = dosis_is - dosis_recomendada_is

    if dosis_is < dosis_recomendada_is:
        estado_dosis_is = "Sub-dosificado"
    else:
        estado_dosis_is = "Sobre-dosificado"

    ahorro_anual_is = diferencia_dosis_is * precio_is * 30 * 12

    if ahorro_anual_is > 0:
        resultado_is = str(
            'Ahorro potencial ' + "%.2f" % ahorro_anual_is) + ' USD/año' + ' al disminuir dosis de antiescala'
    else:
        resultado_is = 'Riesgo de incrustaciones. Aumentar dosis de antiescala'

    # Guardar los resultados de la dosis de antiescala en data frames
    df30.append(dosis_is)
    df31.append(dosis_recomendada_is)
    df32.append(ahorro_anual_is)
    df41.append(estado_dosis_is)

    results_is = pd.DataFrame({'Dosis actual de antiescala [gal/dia]': df30,
                               'Dosis recomendada de antiescala [gal/dia]': df31,
                               'Estado de la dosificacion de antiescala': df41,
                               'Ahorro por optimizacion de antiescala [USD/año]': df32})

    ahorro_total = ahorro_anual_ic + ahorro_anual_is

    corr_button=st.expander('Resultados de corrosion')

    if "corr_state" not in st.session_state:
        st.session_state.corr_state = False

    if corr_button or st.session_state.corr_state:
        st.session_state.corr_state=True
        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric('Velocidad de corrosion en cabeza', str("%.2f" % np.float_(df1[0])) + ' MPY')
            st.metric('Riesgo de corrosion en cabeza', df20[0])

        with col2:
            st.metric('Velocidad de corrosion en fondo', str("%.2f" % np.float_(df2[0])) + ' MPY')
            st.metric('Riesgo de corrosion en fondo', df21[0])

        with col3:
            st.metric('Velocidad de corrosion maxima', str("%.2f" % np.float_(corr_max)) + ' MPY')
            st.metric('Riesgo de corrosion maximo', risk_corr_max)

        st.plotly_chart(fig_corr, use_container_width=True)

        st.write("")

    scale_button=st.expander('Resultados indice de saturacion')

    if "scale_state" not in st.session_state:
        st.session_state.scale_state = False

    if scale_button or st.session_state.scale_state:
        st.session_state.scale_state=True

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric('Indice de saturacion en cabeza', str("%.2f" % np.float_(df3[0])))
            st.metric('Riesgo de incrustaciones en cabeza', df22[0])

        with col2:
            st.metric('Indice de saturacion en fondo', str("%.2f" % np.float_(df4[0])))
            st.metric('Riesgo de incrustaciones en fondo', df23[0])

        with col3:
            st.metric('Indice de saturacion maximo', str("%.2f" % np.float_(scale_max)))
            st.metric('Riesgo de incrustaciones maximo', risk_scale_max)

        st.plotly_chart(fig_sca, use_container_width=True)

    critic_button=st.expander('Calcular criticidad del pozo')

    if "critic_state" not in st.session_state:
        st.session_state.critic_state = False

    if critic_button or st.session_state.critic_state:
        st.session_state.critic_state=True

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric('Produccion de petroleo del pozo', str("%.0f" % np.float_(BOPD)) + ' BOPD')

        with col2:
            st.metric('Indice de criticidad total del pozo', str("%.2f" % np.float_(df9[0])))

        with col3:
            st.metric('Prioridad para tratamiento quimico', nivel_critic)

        st.success('El indice de criticidad se calcula a partir de la produccion de petroleo del pozo, la velocidad de corrosion maxima, y el indice de\
                   saturacion maximo en todo el perfil. Por tanto, es un valor que toma en cuenta el riesgo global de la perdida de produccion en caso de sufrir eventos de\
                   corrosion o incrustaciones a lo largo de toda la tuberia de produccion. El valor del indice de criticidad permite priorizar la dosificacion\
                   de quimicos y el monitoreo a los pozos mas criticos del campo. De esta manera, se optimiza el tiempo y los recursos del personal de\
                   tratamiento quimico.')

    opt_button=st.expander('Optimizar dosis de quimicos')

    if "opt_state" not in st.session_state:
        st.session_state.opt_state = False

    if opt_button or st.session_state.opt_state:
        st.session_state.opt_state=True

        # Dosis anticorrosivo
        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric('Dosis actual anticorrosivo', str("%.0f" % np.float_(dosis_ic)) + ' gal/dia')

        with col2:
            st.metric('Dosis recomendada anticorrosivo', str(dosis_recomendada_ic) + ' gal/dia')

        with col3:
            st.metric('Estado inyeccion anticorrosivo', estado_dosis_ic)

        st.success(resultado_ic)

        # Dosis antiescala
        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric('Dosis actual antiescala', str("%.0f" % np.float_(dosis_is)) + ' gal/dia')

        with col2:
            st.metric('Dosis recomendada antiescala', str(dosis_recomendada_is) + ' gal/dia')

        with col3:
            st.metric('Estado inyeccion antiescala', estado_dosis_is)

        st.success(resultado_is)

        st.markdown('Nota: se asume un precio de 20 USD/gal para el anticorrosivo y el antiescala')

if __name__ == '__main__':
    run()

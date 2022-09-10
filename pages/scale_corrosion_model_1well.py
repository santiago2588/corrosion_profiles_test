import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from funciones import *

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

    add_selectbox = st.sidebar.selectbox(
        "Que deseas predecir?",
        ("Un solo pozo", "Varios pozos"))

    st.title("Cálculo de velocidad de corrosión, índice de saturacion, y criticidad de pozos petroleros")

    if add_selectbox == 'Un solo pozo':

        st.subheader("Predicciones para un solo pozo")

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

        if st.button('Calcular'):
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

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric('Produccion de petroleo del pozo', str("%.0f" % np.float_(BOPD)) + ' BOPD')

            with col2:
                st.metric('Indice de criticidad total del pozo', str("%.2f" % np.float_(df9[0])))

            with col3:
                st.metric('Prioridad para tratamiento quimico', nivel_critic)

            st.success('El indice de criticidad se calcula a partir de la produccion de petroleo del pozo, la velocidad de corrosion promedio, y el indice de\
                       saturacion promedio. Por tanto, es un valor que toma en cuenta el riesgo global de la perdida de produccion en caso de sufrir eventos de\
                       corrosion o incrustaciones a lo largo de toda la tuberia de produccion. El valor del indice de criticidad permite priorizar la dosificacion\
                       de quimicos y el monitoreo a los pozos mas criticos del campo. De esta manera, se optimiza el tiempo y los recursos del personal de\
                       tratamiento quimico.')

        if st.button('Optimizar dosis de quimicos'):
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

        if st.button('Beneficios'):
            st.success(
                "Felicitaciones, has comenzado a desbloquear algunos de los beneficios de la digitalización y el procesamiento de datos con modelos cientificos e inteligencia artificial.")
            st.markdown(
                "En esta ocasion, te hemos presentado los beneficios del modulo de Calculos. Sin embargo, quedan inmensas oportunidades para optimizar la productividad y el desempeño de los pozos petroleros al integrar todos los modulos de ASTRO")
            st.markdown(
                "Si estás listo para emprender este camino y llevar la digitalización al siguiente nivel, contactate con nosotros para guiarte en el proceso y recibir asesoria de nuestros expertos para optimizar tus operaciones.")

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

            local_css("pages/style.css")

    # Codigo para varios pozos

    if add_selectbox == 'Varios pozos':

        st.subheader("Predicciones para múltiples pozos")

        file_upload = st.file_uploader("Subir archivo CSV para las predicciones", type=["csv"])

        if file_upload is not None:
            data = pd.read_csv(file_upload)

            parameters = dict()
            parameters_names = np.array(
                ["BOPD", "BWPD", "MSCF", "Temperature Head", "Temperature Bottom", "Pressure Head", "Pressure Bottom",
                 "CO2 Fraction in Gas", "Alkalinity", "Chlorides", "Sodium", "Magnesium", "Potassium", "Calcium",
                 "Strontium",
                 "Barium", "Sulphates", "Carboxylic acids", "Well Depth", "Well Pipe Diameter",
                 "Corrosion Inhibitor Efficiency",
                 "aux", "dosis_ic", "dosis_is", "precio_ic", "precio_is"])

            for i in data.drop(["Well", "Unit"], axis=1).columns:
                for id_j, j in enumerate(parameters_names):
                    parameters[j] = float(data[i].iloc[int(np.where(data["Well"] == j)[0])])

                # Velocidad de corrosion en cabeza
                nk_temp, corr_ic_temp, corr_risk_temp = calcNorsok(parameters["Temperature Head"],
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

                # Velocidad de corrosion en fondo
                nk_temp1, corr_ic_temp1, corr_risk_temp1 = calcNorsok(parameters["Temperature Bottom"],
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

                # Perfil de la velocidad de corrosion
                temp_array, press_array, depth_array, fy_df, ph_df, nk_df, corr_profile_risk = grahpNorskok(
                    parameters["Temperature Head"],
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

                # Indice de saturacion en cabeza
                calcite_si_temp, solid_temp, scale_risk_temp = calCalcite(parameters["Pressure Head"],
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

                # Indice de saturacion en fondo
                calcite_si_temp1, solid_temp1, scale_risk_temp1 = calCalcite(parameters["Pressure Bottom"],
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

                # Perfil del indice de saturacion
                temp_array, press_array, depth_array, fy, ph1, calcite, ptb1, scale_profile_risk = graphCalcite(
                    parameters["Temperature Head"],
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
                df5.append(parameters['BOPD'])

                # Guardar los resultados del perfil de velocidad de corrosion en un data frame
                df10.append(temp_array)
                df11.append(press_array)
                df12.append(depth_array)
                df13.append(fy_df)
                df14.append(ph_df)
                df15.append(nk_df)
                df25.append(corr_profile_risk)

                # Guardar los resultados del perfil del indice de saturacion en un data frame
                df16.append(fy)
                df17.append(ph1)
                df18.append(calcite)
                df19.append(ptb1)
                df26.append(scale_profile_risk)

                # Calculo de la criticidad

                prod = parameters['BOPD']

                # Velocidad de corrosion promedio
                for list in df15:
                    corr_max = np.max(list)

                # Indice de saturacion promedio
                for list in df18:
                    scale_max = np.max(list)

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
                df6.append(critic_prod)
                df7.append(critic_corr)
                df8.append(critic_si)
                df9.append(critic_tot)
                df24.append(nivel_critic)
                df44.append(corr_max)
                df45.append(scale_max)
                df46.append(risk_corr_max)
                df47.append(risk_scale_max)

                # Calculo de la dosis de quimico recomendada y el ahorro
                BWPD = parameters['BWPD']
                dosis_ic = parameters["dosis_ic"]
                dosis_is = parameters["dosis_is"]
                precio_ic = parameters["precio_ic"]
                precio_is = parameters["precio_is"]

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

                # Guardar los resultados de dosis recomendadas y ahorro en un data frame
                df33.append(dosis_ic)
                df34.append(dosis_recomendada_ic)
                df35.append(resultado_ic)
                df36.append(dosis_is)
                df37.append(dosis_recomendada_is)
                df38.append(resultado_is)
                df42.append(estado_dosis_ic)
                df43.append(estado_dosis_is)

                # Resultados de corrosion

                results_corr = pd.DataFrame(
                    {'Pozo': df0, 'Velocidad de corrosion cabeza [mpy]': df1, 'Riesgo de corrosion cabeza': df20,
                     'Velocidad de corrosion fondo [mpy]': df2, 'Riesgo de corrosion fondo': df21,
                     'Velocidad de corrosion promedio [mpy]': df44, 'Riesgo de corrosion promedio': df46}).set_index(
                    ['Pozo'])

                results_corr_profile = pd.DataFrame({'Pozo': df0, 'Temperatura [F]': df10, 'Presion [psi]': df11,
                                                     'Profundidad [ft]': df12, 'Fugacidad CO2': df13,
                                                     'pH': df14, 'Velocidad de corrosion (mpy)': df15,
                                                     'Riesgo de corrosion': df25}).set_index(['Pozo']).apply(
                    pd.Series.explode).reset_index()

                for i, (Pozo, subdf) in enumerate(results_corr_profile.groupby('Pozo'), 1):
                    locals()[f'well_corr{i}'] = subdf

                # Resultados de escala

                results_scale = pd.DataFrame(
                    {'Pozo': df0, 'Indice de saturacion cabeza': df3, 'Riesgo de incrustaciones cabeza': df22,
                     'Indice de saturacion fondo': df4, 'Riesgo de incrustaciones fondo': df23,
                     'Indice de saturacion promedio [mpy]': df45, 'Riesgo de incrustaciones promedio': df47}).set_index(
                    ['Pozo'])

                results_scale_profile = pd.DataFrame({'Pozo': df0, 'Temperatura [F]': df10, 'Presion [psi]': df11,
                                                      'Profundidad [ft]': df12, 'Fugacidad CO2': df16,
                                                      'pH': df17, 'Indice de saturacion calcita': df18,
                                                      'Solidos [PTB]': df19,
                                                      'Riesgo de incrustaciones': df26}).set_index(['Pozo']).apply(
                    pd.Series.explode).reset_index()

                for i, (Pozo, subdf) in enumerate(results_scale_profile.groupby('Pozo'), 1):
                    locals()[f'well_scale{i}'] = subdf

                # Resultados criticidad

                results_critic = pd.DataFrame(
                    {'Pozo': df0, 'Producción [bopd]': df5, 'Velocidad de corrosion maxima [mpy]': df44,
                     'Indice de saturacion maximo [mpy]': df45,
                     'Criticidad total': df9, 'Prioridad TQ': df24})

                # Resultados optimizacion quimicos

                results_opt_corr = pd.DataFrame({'Pozo': df0, 'Dosis actual de anticorrosivo [gal/dia]': df33,
                                                 'Dosis recomendada de anticorrosivo [gal/dia]': df34,
                                                 'Estado dosificacion de anticorrosivo': df42,
                                                 'Analisis inyeccion de anticorrosivo': df35})

                results_opt_scale = pd.DataFrame({'Pozo': df0, 'Dosis actual de antiescala [gal/dia]': df36,
                                                  'Dosis recomendada de antiescala [gal/dia]': df37,
                                                  'Estado dosificacion de antiescala': df43,
                                                  'Analisis inyeccion de antiescala': df38})

            if st.button('Resultados de corrosion'):

                st.dataframe(results_corr)

                def convert_df(df):
                    return df.to_csv().encode('utf-8')

                csv_corr = convert_df(results_corr)
                st.download_button("📥Press to Download", csv_corr, "file.csv", "text/csv", key='download-csv')

                corr_sliced = [v for k, v in results_corr_profile.groupby('Pozo')]

                for df in corr_sliced:
                    fig_corr = px.line(df, x='Velocidad de corrosion (mpy)', y='Profundidad [ft]',
                                       title='Perfil de velocidad de corrosion',
                                       hover_name='Pozo',
                                       hover_data=['Presion [psi]', 'Temperatura [F]', 'Riesgo de corrosion'])

                    fig_corr.update_traces(mode="markers+lines")
                    fig_corr.update_xaxes(showspikes=True, spikecolor='black')
                    fig_corr.update_yaxes(showspikes=True, spikecolor='black')
                    fig_corr.update_yaxes(autorange="reversed")
                    st.plotly_chart(fig_corr, use_container_width=True)

            if st.button('Resultados indice de saturacion'):

                st.dataframe(results_scale)

                def convert_df(df):
                    return df.to_csv().encode('utf-8')

                csv_scale = convert_df(results_scale)
                st.download_button("📥Press to Download", csv_scale, "file.csv", "text/csv", key='download-csv')

                scale_sliced = [v for k, v in results_scale_profile.groupby('Pozo')]

                for df in scale_sliced:
                    fig_sca = px.line(df, x='Indice de saturacion calcita', y='Profundidad [ft]',
                                      title='Perfil del indice de saturacion',
                                      hover_name='Pozo',
                                      hover_data=['Presion [psi]', 'Temperatura [F]', 'Solidos [PTB]',
                                                  'Riesgo de incrustaciones'])

                    fig_sca.update_traces(mode="markers+lines")
                    fig_sca.update_xaxes(showspikes=True, spikecolor='black')
                    fig_sca.update_yaxes(showspikes=True, spikecolor='black')
                    fig_sca.update_yaxes(autorange="reversed")
                    st.plotly_chart(fig_sca, use_container_width=True)

            if st.button('Criticidad de pozos'):
                st.dataframe(results_critic)

                def convert_df(df):
                    return df.to_csv().encode('utf-8')

                csv_critic = convert_df(results_critic)
                st.download_button("📥Press to Download", csv_critic, "file.csv", "text/csv", key='download-csv')

                fig_crit = px.bar(results_critic, x='Pozo', y='Criticidad total', hover_data=['Prioridad TQ'])
                fig_crit.update_layout(xaxis={'categoryorder': 'total descending'}, title='Criticidad total Pozos')
                st.plotly_chart(fig_crit, use_container_width=True)

                y1 = np.float_(df6)
                y2 = np.float_(df7)
                y3 = np.float_(df8)
                ytot = y1 + y2 + y3

                fig_con = go.Figure()
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y1 / ytot) * 100, name='Produccion'))
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y2 / ytot) * 100, name='Corrosion'))
                fig_con.add_trace(go.Bar(x=results_critic['Pozo'], y=(y3 / ytot) * 100, name='Escala'))
                fig_con.update_layout(barmode='stack', title='Contribuciones (%) a la Criticidad total',
                                      yaxis_title='Contribucion, %')
                st.plotly_chart(fig_con, use_container_width=True)

            if st.button('Optimizar dosis de quimicos'):
                st.dataframe(results_opt_corr)

                def convert_df(df):
                    return df.to_csv().encode('utf-8')

                csv_opt_corr = convert_df(results_opt_corr)
                st.download_button("📥Press to Download", csv_opt_corr, "file.csv", "text/csv", key='download-csv')

                st.dataframe(results_opt_scale)

                def convert_df(df):
                    return df.to_csv().encode('utf-8')

                csv_opt_scale = convert_df(results_opt_scale)
                st.download_button("📥Press to Download", csv_opt_scale, "file.csv", "text/csv", key='download-csv')


# In[20]:


if __name__ == '__main__':
    run()

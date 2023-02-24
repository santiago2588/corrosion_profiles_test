import math

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from typing import List

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
df48=[]
df49=[]


# Calculo de los resultados

# In[19]:
def run():
    st.set_page_config(layout="wide")

    st.title('Plan Estandar')
    st.subheader("Cálculo de velocidad de corrosión, índice de saturacion, y criticidad de pozos petroleros para multiples pozos petroleros")

    st.write('Por favor, sigue los pasos que se presentan a continuacion')
    st.write("1. Descarga la plantilla Excel y llena los datos requeridos de los pozos petroleros, en las unidades correspondientes")
    st.write("2. Carga el archivo Excel lleno")
    st.write("3. Presiona el boton Calcular para ver los resultados: velocidad de corrosion, indice de saturacion y criticidad de los pozos")
    st.write("4. Presiona el boton Optimizar para obtener las dosis recomendadas de los quimicos, junto con el ahorro generado")

    filepath='Databases/plantilla_pozos.csv'
    with open(filepath, 'rb') as excel_template:
        st.download_button(label = 'Descargar plantilla Excel', data = excel_template, file_name = 'plantilla_pozos.csv', mime = 'text/csv')

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

            # Velocidad de corrosion maximo
            for list in df15:
                corr_max = np.max(list)

            # Indice de saturacion maximo
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

            #Funcion para redondear a 0.5 la dosis recomendada. No la utilizo para tener un rango de seguridad en la inyeccion, porque utilizo math.ceil que redondea el numero al entero superior
            def round_to_half(n):
                return round(n*2)/2

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
                resultado_ic = 'Ahorro potencial al disminuir dosis de anticorrosivo'
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
                resultado_is = 'Ahorro potencial al disminuir dosis de antiescala'
            else:
                resultado_is = 'Riesgo de incrustaciones. Aumentar dosis de antiescala'

            # Guardar los resultados de dosis recomendadas y ahorro en un data frame
            df33.append(dosis_ic)
            df34.append(dosis_recomendada_ic)
            df35.append(resultado_ic)
            df48.append(ahorro_anual_ic)
            df36.append(dosis_is)
            df37.append(dosis_recomendada_is)
            df38.append(resultado_is)
            df42.append(estado_dosis_ic)
            df43.append(estado_dosis_is)
            df49.append(ahorro_anual_is)

            # Resultados de corrosion

            results_corr = pd.DataFrame(
                {'Pozo': df0, 'Velocidad de corrosion cabeza [mpy]': df1, 'Riesgo de corrosion cabeza': df20,
                 'Velocidad de corrosion fondo [mpy]': df2, 'Riesgo de corrosion fondo': df21,
                 'Velocidad de corrosion maximo [mpy]': df44, 'Riesgo de corrosion maximo': df46}).set_index(
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
                 'Indice de saturacion maximo [mpy]': df45, 'Riesgo de incrustaciones maximo': df47}).set_index(
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
                                             'Analisis inyeccion de anticorrosivo': df35,
                                             'Ahorro potencial [USD/año]':df48})

            results_opt_scale = pd.DataFrame({'Pozo': df0, 'Dosis actual de antiescala [gal/dia]': df36,
                                              'Dosis recomendada de antiescala [gal/dia]': df37,
                                              'Estado dosificacion de antiescala': df43,
                                              'Analisis inyeccion de antiescala': df38,
                                              'Ahorro potencial [USD/año]':df49})

        st.write("## Resultados")

        with st.expander('Resultados de corrosion'):

            st.dataframe(results_corr)

            def convert_df(df):
                return df.to_csv().encode('utf-8')

            csv_corr = convert_df(results_corr)
            st.download_button("📥Press to Download", csv_corr, "file.csv", "text/csv", key='download-csv')

            st.write("Perfiles velocidad de corrosion")

            corr_sliced = [v for k, v in results_corr_profile.groupby('Pozo')]
            for df in corr_sliced:
                fig_corr = px.line(df, x='Velocidad de corrosion (mpy)', y='Profundidad [ft]',
                                   hover_name='Pozo',
                                   hover_data=['Presion [psi]', 'Temperatura [F]', 'Riesgo de corrosion'])

                fig_corr.update_traces(mode="markers+lines")
                fig_corr.update_xaxes(showspikes=True, spikecolor='black')
                fig_corr.update_yaxes(showspikes=True, spikecolor='black')
                fig_corr.update_yaxes(autorange="reversed")
                st.plotly_chart(fig_corr, use_container_width=True)

        with st.expander('Resultados indice de saturacion'):

            st.dataframe(results_scale)

            def convert_df(df):
                return df.to_csv().encode('utf-8')

            csv_scale = convert_df(results_scale)
            st.download_button("📥Press to Download", csv_scale, "file.csv", "text/csv", key='download-csv1')

            st.write("Perfiles indice de saturacion")

            def plot_figures(figures: List[go.Figure], tab_names: List[str]):
                tabs=st.tabs(tab_names)
                for i, (fig,tab_name) in enumerate(zip(figures, tab_names)):
                    with tabs[i]:
                        st.plotly_chart(fig)

            scale_sliced = [v for k, v in results_scale_profile.groupby('Pozo')]
            st.write(scale_sliced)
            #st.write(scale_sliced.type())

            dfx=scale_sliced[1]
            st.write(dfx)

            tab_names=df0

            # for df in scale_sliced:
            #     fig_sca = px.line(df, x='Indice de saturacion calcita', y='Profundidad [ft]',
            #                       hover_name='Pozo',
            #                       hover_data=['Presion [psi]', 'Temperatura [F]', 'Solidos [PTB]',
            #                                   'Riesgo de incrustaciones'])
            #
            #     fig_sca.update_traces(mode="markers+lines")
            #     fig_sca.update_xaxes(showspikes=True, spikecolor='black')
            #     fig_sca.update_yaxes(showspikes=True, spikecolor='black')
            #     fig_sca.update_yaxes(autorange="reversed")
            #
            #     #tab=st.tabs(df0)
            #     #with tab:
            #     st.plotly_chart(fig_sca, use_container_width=True)
            #     #plot_figures([fig_sca],df0)

            for i, df in enumerate(scale_sliced):
                with st.tab(tab_names[i]):
                    fig_sca = px.line(df, x='Indice de saturacion calcita', y='Profundidad [ft]',
                                           hover_name='Pozo',
                                           hover_data=['Presion [psi]', 'Temperatura [F]', 'Solidos [PTB]',
                                                       'Riesgo de incrustaciones'])

                    fig_sca.update_traces(mode="markers+lines")
                    fig_sca.update_xaxes(showspikes=True, spikecolor='black')
                    fig_sca.update_yaxes(showspikes=True, spikecolor='black')
                    fig_sca.update_yaxes(autorange="reversed")
                    st.plotly_chart(fig_sca)


        with st.expander('Criticidad de pozos'):

            st.dataframe(results_critic)

            def convert_df(df):
                return df.to_csv().encode('utf-8')

            csv_critic = convert_df(results_critic)
            st.download_button("📥Press to Download", csv_critic, "file.csv", "text/csv", key='download-csv2')

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

        with st.expander('Optimizar dosis de quimicos'):

            def convert_df(df):
                return df.to_csv().encode('utf-8')

            st.write("Optimizacion Anticorrosivo")
            st.dataframe(results_opt_corr)
            csv_opt_corr = convert_df(results_opt_corr)
            st.download_button("📥Press to Download", csv_opt_corr, "file.csv", "text/csv", key='download-csv3')
            st.metric(label='Ahorro total Anticorrosivo [USD/año]',value=sum(results_opt_corr['Ahorro potencial [USD/año]']))

            st.write("Optimizacion Antiescala")
            st.dataframe(results_opt_scale)
            csv_opt_scale = convert_df(results_opt_scale)
            st.download_button("📥Press to Download", csv_opt_scale, "file.csv", "text/csv", key='download-csv4')
            st.metric(label='Ahorro total Antiescala [USD/año]',value=sum(results_opt_scale['Ahorro potencial [USD/año]']))

# In[20]:


if __name__ == '__main__':
    run()

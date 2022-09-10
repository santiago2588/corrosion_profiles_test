import streamlit as st
from PIL import Image

st.set_page_config(layout="wide")

st.markdown("# Informacion")

st.write("""
        ### Bienvenido al demo de ASTRO: la herramienta digital para la gestión del tratamiento químico en la industria petrolera
        """)

with st.expander('Para que sirve ASTRO❓'):
    
    st.markdown("""
    ### ASTRO genera los siguientes beneficios:
    - Reduce los costos asociados al tratamiento químico y las pérdidas de producción por eventos de corrosión e incrustaciones
    - Identifica los pozos críticos en tus operaciones con nuestra metodología para priorizar y optimizar el tratamiento químico
    - Disminuye el tiempo para el procesamiento de datos, genera reportes automáticos y libera tiempo valioso para optimizar la rentabilidad y seguridad de las operaciones""")

    st.write("")

    st.markdown("ASTRO es una solucion modular que se adapta a tus necesidades y genera valor en tus operaciones. ASTRO incluye los siguientes modulos:")
    
    image1 = Image.open('Resources/modulos astro.png')
    st.image(image1)
      
    st.markdown('En esta demostracion, te presentamos el Modulo de Calculos, que permite calcular la velocidad de corrosion, indice de saturacion y criticidad de los pozos petroleros.')
                 
with st.expander("Para que sirve el modulo de calculos de ASTRO❓"):

    st.markdown("""
    ### Este modulo permite predecir la velocidad de corrosión, índice de saturacion, y criticidad de los pozos de petróleo en función de varios parametros de produccion y laboratorio. El modulo de calculos tambien permite optimizar las dosis de los quimicos anticorrosivo y antiescala en funcion de los resultados
    obtenidos en el modelo (riesgo de corrosion, incrustaciones y criticidad). 
    - Si deseas obtener las predicciones para 1 pozo, selecciona la hoja "Calculos 1 pozo". Luego, ajusta el valor de cada parametro y presiona el boton Calcular para ver los resultados. Despues, presiona el boton Optimizar y obtendras las dosis recomendadas de los quimicos, junto con el ahorro generado.
    - Si deseas obtener las predicciones para multiples pozos, selecciona la hoja "Calculos varios pozos". Carga el archivo CSV con los parametros de produccion y laboratorio de los pozos. Despues, da click en los botones que se indican para obtener los resultados. 
    """)

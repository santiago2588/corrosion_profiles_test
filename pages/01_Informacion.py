import streamlit as st
from PIL import Image

st.set_page_config(layout="wide")

st.markdown("# Informacion")

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
    ### El modulo de Calculos de ASTRO te sirve para:
    - Predecir la velocidad de corrosión, índice de saturacion, y criticidad de los pozos de petróleo en función de varios parametros de produccion y laboratorio. 
    - Optimizar las dosis de los quimicos anticorrosivo y antiescala en funcion de los resultados obtenidos en el modelo (riesgo de corrosion, incrustaciones y criticidad).
    - Generar ahorros por optimizacion de los quimicos e incremento de la produccion de los pozos debido a un tratamiento quimico adecuado. 
    """)

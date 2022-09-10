# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 18:19:54 2022

@author: mjkipsz2
"""

import streamlit as st

st.set_page_config(layout="wide",page_title="ASTRO",page_icon="🌿")

from PIL import Image
image = Image.open('Resources/logo_Pungo.png')
st.image(image)

st.write("# DIGITALIZACIÓN QUE GENERA IMPACTO")
st.write("### Soluciones digitales para incrementar tus ingresos, disminuir tus costos operativos, reducir riesgos en tus operaciones y mejorar tu desempeño ambiental, desarrolladas a la medida de tus necesidades.")

st.write("""
        ### Bienvenido al demo de ASTRO: la herramienta digital para la gestión del tratamiento químico en la industria petrolera:
        - Reduce los costos asociados al tratamiento químico y las pérdidas de producción por eventos de corrosión e incrustaciones
        - Identifica los pozos críticos en tus operaciones con nuestra metodología para priorizar y optimizar el tratamiento químico
        - Disminuye el tiempo para el procesamiento de datos, genera reportes automáticos y libera tiempo valioso para optimizar la rentabilidad y seguridad de las operaciones
        """)


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
st.write("### Bienvenido al demo de ASTRO: la herramienta digital para la gestion del tratamiento quimico en la industria petrolera, desarrollada a la medida de tus necesidades.")




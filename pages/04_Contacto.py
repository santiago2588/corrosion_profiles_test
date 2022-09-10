import streamlit as st

st.markdown("# Felicitaciones, has comenzado a desbloquear algunos de los beneficios de la digitalización y el procesamiento de datos con modelos cientificos e inteligencia artificial.")

st.write(""" En esta ocasion, te hemos presentado los beneficios del modulo de Calculos. Sin embargo, quedan inmensas oportunidades para optimizar la productividad y el desempeño de los pozos petroleros al integrar todos los modulos de ASTRO.
    
    Si estás listo para emprender este camino y llevar la digitalización al siguiente nivel, contactate con nosotros para guiarte en el proceso y recibir asesoria de nuestros expertos para optimizar tus operaciones.
""")

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

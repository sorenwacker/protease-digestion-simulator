import streamlit as st

from digest_simulator.DigestionSimulator import DigestionSimulator  # assuming the class is defined in digestion_simulator.py
from digest_simulator.proteases import available_proteases

st.set_page_config(layout="wide")

# Use the write function for basic output
st.title('Protein Digestion Simulator')

hemoglobin_beta_sequence = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'

# Provide the sequence as an input with a default value
sequence = st.text_input('Enter protein sequence', value=hemoglobin_beta_sequence)

# Add a multiselect field to choose proteases
selected_proteases = st.multiselect('Select proteases', list(available_proteases.keys()))

# Initialize the DigestionSimulator with the entered sequence and selected proteases
min_peptide_length = st.slider('Minimum peptide length:', min_value=0, max_value=10, value=3, step=1)


if sequence and selected_proteases:

    st.header('Digestion prediction')    
    # Map the names of the selected proteases to their instances
    selected_protease_instances = [available_proteases[protease]() for protease in selected_proteases]

    simulator = DigestionSimulator(sequence, selected_protease_instances, min_peptide_length=min_peptide_length)

    # Display the peptide tree
    st.subheader('Peptide Tree')
    
    st.text(simulator.draw_tree().replace(' ', '-'))

    # Extract and display unique peptide sequences
    st.subheader('Predicted (unique) peptide sequences.')
    simulator.extract_unique_peptide_sequences()
    st.text('\n'.join(sorted(simulator.unique_peptide_sequences, key=len)))

    st.header('Predict proteases')
    # Input field for user to add their own sequences
    user_sequences = st.text_area('Enter your own experimentally observed sequences (each on a new line).')

    # Identify and display known proteases
    st.subheader('Predicted proteases from observed sequences.')
    if user_sequences:
        # convert the user input into a list of sequences
        user_sequences = user_sequences.split('\n')
        df = simulator.predict_proteases(user_sequences)
        df['Probabilty [%]'] = df['Probability'] * 100
        df = df.sort_values('Probabilty [%]', ascending=False).reset_index(drop=True)
        df.drop('Probability', axis=1, inplace=True)
        st.markdown(df.round(0).to_markdown())

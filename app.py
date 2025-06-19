import streamlit as st
import pandas as pd
import numpy as np
from datetime import datetime
import os
from dotenv import load_dotenv
from pymongo import MongoClient
import re
import requests
import json
import plotly.express as px
import plotly.graph_objects as go
from streamlit_autorefresh import st_autorefresh
import py3Dmol
import streamlit.components.v1 as components

# Load environment variables
load_dotenv()

# Constants
API_URL = "http://localhost:8000"  # FastAPI backend URL

# Page config
st.set_page_config(
    page_title="Mutalyze - Gene Mutation Impact Predictor",
    page_icon="🧬",
    layout="wide"
)

# Initialize session state
if 'user' not in st.session_state:
    st.session_state.user = None
if 'last_prediction' not in st.session_state:
    st.session_state.last_prediction = None

# Authentication functions
def login_user(email: str, password: str):
    # In production, use proper authentication service
    # This is a mock implementation
    if email == "demo@mutalyze.com" and password == "demo123":
        st.session_state.user = {"email": email, "name": "Demo User"}
        return True
    return False

def logout_user():
    st.session_state.user = None

# Sidebar
with st.sidebar:
    st.title("🧬 Mutalyze")
    
    # Authentication
    if st.session_state.user is None:
        st.header("Login")
        email = st.text_input("Email")
        password = st.text_input("Password", type="password")
        if st.button("Login"):
            if login_user(email, password):
                st.success("Login successful!")
                st.experimental_rerun()
            else:
                st.error("Invalid credentials")
    else:
        st.write(f"Welcome, {st.session_state.user['name']}")
        if st.button("Logout"):
            logout_user()
            st.experimental_rerun()
    
    st.header("About")
    st.write("""
    Mutalyze helps predict the impact of gene mutations using AI models.
    Enter a gene name and mutation to get started.
    """)
    
    st.header("Example")
    st.code("Gene: BRCA1\nMutation: A123T")
    
    # Batch upload option
    st.header("Batch Analysis")
    uploaded_file = st.file_uploader("Upload CSV file", type=['csv'])
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        if st.button("Analyze Batch"):
            with st.spinner("Analyzing mutations..."):
                # Prepare batch request
                mutations = [{"gene_name": row['gene'], "mutation": row['mutation']} 
                           for _, row in df.iterrows()]
                
                # Call API
                response = requests.post(f"{API_URL}/batch-predict", 
                                      json={"mutations": mutations})
                if response.status_code == 200:
                    st.session_state.batch_results = response.json()
                    st.success("Batch analysis complete!")
                else:
                    st.error("Error in batch analysis")

# Main content
if st.session_state.user is None:
    st.warning("Please login to use Mutalyze")
    st.stop()

# Main content in a single column
st.header("Input")
gene_name = st.text_input("Gene Name", placeholder="e.g., BRCA1")
mutation = st.text_input("Mutation", placeholder="e.g., A123T")

# Validate mutation format
def validate_mutation(mutation):
    pattern = r'^[A-Z]\d+[A-Z]$'
    return bool(re.match(pattern, mutation))

if st.button("Predict Impact"):
    if not gene_name or not mutation:
        st.error("Please enter both gene name and mutation")
    elif not validate_mutation(mutation):
        st.error("Invalid mutation format. Please use format like 'A123T'")
    else:
        with st.spinner("Analyzing mutation..."):
            # Call API
            response = requests.post(f"{API_URL}/predict",
                                  json={"gene_name": gene_name, "mutation": mutation})
            if response.status_code == 200:
                st.session_state.last_prediction = response.json()
                st.success("Prediction complete!")
            else:
                st.error(f"Error in prediction: {response.status_code}") # Added status code for debug

# Results section (conditional)
if st.session_state.last_prediction:
    st.header("Results")
    result = st.session_state.last_prediction

    # Display main metrics
    col_metric1, col_metric2 = st.columns(2) # Keep metrics in columns for layout
    with col_metric1:
        st.metric("Predicted Impact", result['impact'])
    with col_metric2:
        st.metric("Confidence Score", f"{result['confidence']:.2%}")

    # Display additional information
    st.subheader("Additional Information")
    st.write(f"""
    - **Protein Name**: {result['protein_name']}
    - **Protein Domain**: {result['protein_domain']}
    - **Conservation Score**: {result['conservation_score']:.2f}
    - **Clinical Significance**: {result['clinical_significance']}
    - **UniProt ID**: {result['uniprot_id']}
    """)

    # Add visualization
    st.subheader("Conservation Score Visualization")
    fig = px.bar(
        x=['Conservation Score'],
        y=[result['conservation_score']],
        range_y=[0, 1],
        title="Mutation Conservation Score"
    )
    st.plotly_chart(fig)

    # Display protein structure if available
    if result['structure_url']:
        st.subheader("Protein Structure")
        try:
            pdb_url = result['structure_url']

            # Create a unique ID for the viewer element
            viewer_id = f"pdb_viewer_{result.get('uniprot_id', 'unknown')}"

            # Fetch PDB data
            pdb_response = requests.get(pdb_url)
            if pdb_response.status_code == 200:
                pdb_data = pdb_response.text

                # Create HTML for 3Dmol viewer with explicit initialization
                html_string = f"""
                <div style=\"height: 400px; width: 100%;\" id=\"{viewer_id}\"></div>
                <script src=\"https://3Dmol.org/build/3Dmol-min.js\"></script>
                <script>
                    console.log('Loading 3Dmol viewer for {viewer_id}');
                    var element = document.getElementById('{viewer_id}');
                    if (element) {{
                        console.log('Element {viewer_id} found. Creating viewer.');
                        var config = {{ backgroundColor: 'white' }};
                        var viewer = 3Dmol.createViewer(element, config);
                        console.log('Viewer created. Adding model.');
                        viewer.addModel(`{pdb_data}`, "pdb");
                        viewer.setStyle({{'cartoon': {{'color': 'spectrum'}}}});
                        viewer.zoomTo();
                        viewer.render();
                        console.log('Viewer rendered.');
                    }} else {{
                        console.log('Element {viewer_id} not found.');
                    }}
                </script>
                """

                # Embed the HTML component
                components.html(html_string, height=400, scrolling=False)

            else:
                st.error(f"Error fetching PDB data from {pdb_url}. Status code: {pdb_response.status_code}")

        except requests.exceptions.RequestException as e:
            st.error(f"Error fetching protein structure: {e}")
        except Exception as e:
            st.error(f"Error rendering protein structure: {e}")

# Batch Results Section
if 'batch_results' in st.session_state:
    st.header("Batch Analysis Results")
    results_df = pd.DataFrame(st.session_state.batch_results)
    st.dataframe(results_df)
    
    # Batch visualization
    st.subheader("Impact Distribution")
    impact_counts = results_df['impact'].value_counts()
    fig = px.pie(
        values=impact_counts.values,
        names=impact_counts.index,
        title="Distribution of Mutation Impacts"
    )
    st.plotly_chart(fig)

# History section
st.header("Prediction History")
try:
    response = requests.get(f"{API_URL}/history")
    if response.status_code == 200:
        history = response.json()
        if history:
            df = pd.DataFrame(history)
            df['timestamp'] = pd.to_datetime(df['timestamp']).dt.strftime('%Y-%m-%d %H:%M:%S')
            st.dataframe(df[['gene', 'mutation', 'impact', 'confidence', 'timestamp']])
        else:
            st.info("No prediction history available")
    else:
        st.warning("Could not fetch prediction history")
except Exception as e:
    st.warning(f"Error connecting to API: {e}")

# Footer
st.markdown("---")

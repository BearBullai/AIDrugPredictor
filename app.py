import streamlit as st
import pandas as pd
import numpy as np
import os
import sys
import logging
from pathlib import Path
from typing import Optional, Dict, Any
import plotly.express as px

# Import local modules
from utils import TargetProcessor, setup_directories, save_analysis_results
from config import config, ensure_dir

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set page config
st.set_page_config(
    page_title="AI Drug Target Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Set up directories
dirs = setup_directories(Path(__file__).parent)
DATA_DIR = Path(dirs['data'])
TARGETS_DIR = Path(dirs['targets'])
STRUCTURES_DIR = Path(dirs['structures'])
POCKETS_DIR = Path(dirs['pockets'])
LIGANDS_DIR = Path(dirs['ligands'])
NANODELIVERY_DIR = Path(dirs['nanodelivery'])
NOTEBOOKS_DIR = Path(dirs['notebooks'])

# Initialize TargetProcessor
target_processor = TargetProcessor(DATA_DIR)

# Initialize session state
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'df' not in st.session_state:
    st.session_state.df = None
if 'filtered_df' not in st.session_state:
    st.session_state.filtered_df = None
if 'processor' not in st.session_state:
    st.session_state.processor = target_processor

# Handle DrugBank viewer route
def handle_drugbank_viewer():
    """Check if the current page is a DrugBank viewer page."""
    if 'drugbank' in st.query_params:
        drug_id = st.query_params['drugbank']
        st.session_state.current_page = f"DrugBank: {drug_id}"
        st.switch_page("pages/7_DrugBank_Viewer.py")
        return True
    return False

# Main app
def main():
    # Check for DrugBank viewer route first
    if handle_drugbank_viewer():
        return
        
    st.sidebar.title("🧬 Navigation")
    
    # Define menu items
    menu_items = [
        "🏠 Home",
        "📤 Data Upload",
        "🔍 Target Analysis",
        "🧬 Structure Prediction",
        "🕳️ Pocket Detection",
        "⚗️ Ligand Docking",
        "📦 Nanocarrier Design",
        "📊 Reports"
    ]
    
    # Add logo and title to sidebar
    st.sidebar.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=80)
    st.sidebar.markdown("# AI Drug Target Platform")
    st.sidebar.markdown("---")
    
    # Navigation
    choice = st.sidebar.radio("Go to", menu_items)
    
    # Show selected page
    st.title(choice.split(" ", 1)[1] if " " in choice else choice)
    st.markdown("---")
    
    # Page routing
    if "Home" in choice:
        show_home()
    elif "Data Upload" in choice:
        show_data_upload()
    elif "Target Analysis" in choice:
        show_target_analysis()
    elif "Structure Prediction" in choice:
        show_structure_prediction()
    elif "Pocket Detection" in choice:
        show_pocket_detection()
    elif "Ligand Docking" in choice:
        show_ligand_docking()
    elif "Nanocarrier Design" in choice:
        show_nanocarrier_design()
    elif "Reports" in choice:
        show_reports()
        
    # Add footer
    st.sidebar.markdown("---")
    st.sidebar.info(
        "ℹ️ AI Drug Target Platform\n"
        "Version 1.0.0\n\n"
        "For support, contact your system administrator."
    )

def show_home():
    st.markdown("""
    ## 🚀 Welcome to the AI Drug Target Platform
    
    A comprehensive platform for identifying, prioritizing, and validating drug targets with AI-powered analysis.
    """)
    
    # Create a grid of feature cards
    col1, col2 = st.columns(2)
    
    with col1:
        with st.container(border=True):
            st.markdown("### 📊 Data Analysis")
            st.markdown("""
            - Upload and preprocess differential gene expression data
            - Interactive visualization of gene expression profiles
            - Statistical analysis and quality control
            """)
            
        st.write("")
        
        with st.container(border=True):
            st.markdown("### 🧬 Structure Prediction")
            st.markdown("""
            - Predict 3D protein structures using ColabFold/AF2
            - Model quality assessment
            - Structure visualization and analysis
            """)
            
        st.write("")
        
        with st.container(border=True):
            st.markdown("### ⚗️ Ligand Docking")
            st.markdown("""
            - Virtual screening of compound libraries
            - Molecular docking and pose prediction
            - Binding affinity estimation
            """)
    
    with col2:
        with st.container(border=True):
            st.markdown("### 🎯 Target Prioritization")
            st.markdown("""
            - Score and rank potential drug targets
            - Integration of multi-omics data
            - Pathway and network analysis
            """)
            
        st.write("")
        
        with st.container(border=True):
            st.markdown("### 🕳️ Binding Pocket Detection")
            st.markdown("""
            - Identify potential drug binding sites
            - Pocket druggability assessment
            - Conservation analysis
            """)
            
        st.write("")
        
        with st.container(border=True):
            st.markdown("### 📦 Nanocarrier Design")
            st.markdown("""
            - Design targeted drug delivery systems
            - Carrier selection and optimization
            - Formulation recommendations
            """)
    
    # Getting Started Section
    st.markdown("---")
    st.markdown("## 🚀 Getting Started")
    
    steps = [
        "1. **Upload Data**: Go to the 'Data Upload' section and upload your differential gene expression data.",
        "2. **Analyze Targets**: Use the 'Target Analysis' section to filter and prioritize potential drug targets.",
        "3. **Predict Structures**: In 'Structure Prediction', generate 3D models of your target proteins.",
        "4. **Find Binding Pockets**: Use 'Pocket Detection' to identify potential drug binding sites.",
        "5. **Dock Ligands**: In 'Ligand Docking', screen compounds against your target proteins.",
        "6. **Design Nanocarriers**: Use the 'Nanocarrier Design' section to create targeted delivery systems.",
        "7. **Generate Reports**: Finally, create comprehensive reports of your findings in the 'Reports' section."
    ]
    
    for step in steps:
        st.markdown(step)
    
    # Add a quick start button
    if st.button("📤 Start by Uploading Data", type="primary"):
        st.switch_page("pages/1_Data_Upload.py")

def show_data_upload():
    st.header("📤 Upload Your Data")
    
    # File uploader
    uploaded_file = st.file_uploader(
        "Upload your DEG list (CSV or Excel)", 
        type=config.get('data.allowed_extensions', ['.csv', '.xlsx', '.xls'])
    )
    
    if uploaded_file is not None:
        try:
            # Process the uploaded file
            df = st.session_state.processor.load_data(uploaded_file)
            df = st.session_state.processor.preprocess_data(df)
            
            # Update session state
            st.session_state.df = df
            st.session_state.data_loaded = True
            st.session_state.filtered_df = None
            
            # Show success message
            st.success("✅ Data loaded and preprocessed successfully!")
            
            # Show data preview
            with st.expander("📊 Data Preview", expanded=True):
                st.dataframe(df.head())
                
                # Basic stats
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Number of Genes", len(df))
                with col2:
                    st.metric("Mean |logFC|", f"{abs(df['logfc']).mean():.2f}")
                with col3:
                    st.metric("Mean p-value", f"{df['adj.p.val'].mean():.3f}")
            
            # Save the uploaded file
            save_path = DATA_DIR / "deg_list.csv"
            df.to_csv(save_path, index=False)
            
            # Show file info
            st.info(f"💾 Data saved to: `{save_path}`")
            
            # Enable the analysis tab
            st.session_state.show_analysis = True
            
        except Exception as e:
            st.error(f"❌ Error processing file: {str(e)}")
            logger.exception("Error in show_data_upload")
    else:
        st.info("ℹ️ Please upload a file to get started")
        
        # Check if data already exists
        default_data = DATA_DIR / "deg_list.csv"
        if default_data.exists():
            if st.button("📂 Load existing data"):
                try:
                    df = pd.read_csv(default_data)
                    df = st.session_state.processor.preprocess_data(df)
                    
                    st.session_state.df = df
                    st.session_state.data_loaded = True
                    st.session_state.filtered_df = None
                    st.session_state.show_analysis = True
                    
                    st.success("✅ Existing data loaded successfully!")
                    st.dataframe(df.head())
                except Exception as e:
                    st.error(f"❌ Error loading existing data: {str(e)}")
                    logger.exception("Error loading existing data")

def show_target_analysis():
    st.header("🔍 Target Analysis and Prioritization")
    
    if not st.session_state.data_loaded:
        st.warning("⚠️ Please upload data in the 'Data Upload' section first.")
        return
    
    df = st.session_state.df
    
    # Get default values from config
    default_min_logfc = config.get('analysis.default_min_logfc', 1.0)
    default_max_pval = config.get('analysis.default_max_pval', 0.05)
    top_n = config.get('analysis.top_n_targets', 50)
    
    # Create two columns for filters and visualization
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("🔧 Filter Settings")
        
        # Get min/max values for sliders
        logfc_min, logfc_max = float(df['logfc'].min()), float(df['logfc'].max())
        pval_min, pval_max = float(df['adj.p.val'].min()), float(df['adj.p.val'].max())
        
        # Ensure we have valid ranges
        if abs(logfc_max - logfc_min) < 0.1:
            logfc_min, logfc_max = -5.0, 5.0
        
        # Create sliders
        logfc_threshold = st.slider(
            "Log2 Fold Change Threshold", 
            min_value=0.0, 
            max_value=max(5.0, abs(logfc_max)), 
            value=default_min_logfc,
            step=0.1,
            help="Minimum absolute log2 fold change for significant genes"
        )
        
        pval_threshold = st.slider(
            "Adjusted p-value Threshold", 
            min_value=0.0, 
            max_value=0.1, 
            value=default_max_pval,
            step=0.001,
            format="%.3f",
            help="Maximum adjusted p-value for significant genes"
        )
        
        # Filter button
        if st.button("Apply Filters", type="primary"):
            try:
                filtered_df = st.session_state.processor.filter_targets(
                    df, 
                    min_logfc=logfc_threshold, 
                    max_pval=pval_threshold
                )
                st.session_state.filtered_df = filtered_df
                st.session_state.filter_params = {
                    'logfc_threshold': logfc_threshold,
                    'pval_threshold': pval_threshold
                }
                st.rerun()
            except Exception as e:
                st.error(f"❌ Error filtering data: {str(e)}")
                logger.exception("Error in target analysis")
    
    with col2:
        st.subheader("📊 Volcano Plot")
        
        # Show volcano plot
        try:
            fig = st.session_state.processor.create_volcano_plot(
                df,
                logfc_threshold=st.session_state.get('filter_params', {}).get('logfc_threshold', default_min_logfc),
                pval_threshold=st.session_state.get('filter_params', {}).get('pval_threshold', default_max_pval)
            )
            if fig:
                st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.error(f"❌ Error generating volcano plot: {str(e)}")
    
    # Show filtered results if available
    if st.session_state.filtered_df is not None and not st.session_state.filtered_df.empty:
        filtered_df = st.session_state.filtered_df
        
        st.subheader(f"🎯 Filtered Targets (n={len(filtered_df)})")
        
        # Show top targets
        st.dataframe(
            filtered_df.head(top_n).style.background_gradient(
                subset=['score'], 
                cmap='YlOrRd'
            ),
            height=400,
            use_container_width=True
        )
        
        # Download button
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Results",
            data=csv,
            file_name="filtered_targets.csv",
            mime="text/csv"
        )
        
        # Show some statistics
        st.subheader("📈 Statistics")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Targets", len(filtered_df))
        with col2:
            st.metric("Mean Score", f"{filtered_df['score'].mean():.2f}")
        with col3:
            st.metric("Top Score", f"{filtered_df['score'].max():.2f}")
    else:
        st.info("ℹ️ Apply filters to see the filtered results")

# Run the app
if __name__ == "__main__":
    main()

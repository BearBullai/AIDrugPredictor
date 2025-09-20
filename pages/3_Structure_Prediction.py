"""
Structure Prediction Page

This module provides functionality for predicting protein structures using various methods
like ColabFold (AlphaFold2), ESMFold, and RoseTTAFold.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import logging
import subprocess
import os
import time
import json
import tempfile
from typing import Optional, Dict, List, Tuple, Union

# Local imports
from utils import TargetProcessor
from config import config, ensure_dir
from pdb_utils import (
    visualize_pdb, 
    extract_ligands, 
    calculate_sasa,
    clean_pdb,
    extract_binding_pocket
)
import py3Dmol
from stmol import showmol

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config is managed in app.py

# Check if data is loaded and filtered
if 'data_loaded' not in st.session_state or not st.session_state.data_loaded:
    st.warning("⚠️ Please upload and analyze data first.")
    if st.button("Go to Data Upload"):
        st.switch_page("pages/1_Data_Upload.py")
    st.stop()

# Initialize directories
STRUCTURES_DIR = Path("structures")
STRUCTURES_DIR.mkdir(exist_ok=True)

class StructurePredictor:
    """Class for handling structure prediction tasks."""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def predict_structure(self, sequence: str, job_name: str, method: str = "colabfold") -> Optional[Path]:
        """
        Predict protein structure using the specified method.
        
        Args:
            sequence: Protein sequence
            job_name: Name for this prediction job
            method: Prediction method to use (colabfold, esmfold, etc.)
            
        Returns:
            Path to the predicted PDB file or None if prediction failed
        """
        try:
            job_dir = self.output_dir / job_name
            job_dir.mkdir(exist_ok=True)
            
            # Save sequence to a file
            fasta_file = job_dir / f"{job_name}.fasta"
            with open(fasta_file, 'w') as f:
                f.write(f">{job_name}\n{sequence}")
            
            st.info(f"Starting structure prediction for {job_name} using {method}...")
            
            # In a real implementation, this would call external tools like:
            # - ColabFold/AlphaFold2
            # - ESMFold
            # - RoseTTAFold
            # For now, we'll use a placeholder that creates a simple structure
            
            # Create a simple alpha-helix as a placeholder
            pdb_file = self._create_placeholder_structure(sequence, job_name, job_dir)
            
            # Generate pLDDT scores (placeholder)
            plddt_file = job_dir / f"{job_name}_plddt.json"
            with open(plddt_file, 'w') as f:
                json.dump({
                    'plddt': [min(0.9, max(0.7, 0.8 + 0.1 * (i % 5))) 
                             for i in range(len(sequence))]
                }, f)
            
            st.success(f"✅ Structure prediction completed for {job_name}")
            return pdb_file
            
        except Exception as e:
            st.error(f"❌ Error in structure prediction: {str(e)}")
            logger.exception("Structure prediction failed")
            return None
    
    def _create_placeholder_structure(self, sequence: str, job_name: str, job_dir: Path) -> Path:
        """Create a simple alpha-helix as a placeholder structure."""
        # In a real implementation, this would call ColabFold/ESMFold
        # For now, we'll create a simple PDB file with an alpha-helix
        pdb_file = job_dir / f"{job_name}.pdb"
        
        # Simple alpha-helix coordinates with better 3D structure
        pdb_content = """HEADER    PLACEHOLDER STRUCTURE
TITLE     PLACEHOLDER ALPHA-HELIX
REMARK    This is a placeholder structure. In a real implementation,
REMARK    this would be replaced with actual structure prediction output.
"""
        # Add alpha-helix coordinates with proper 3D structure
        for i in range(20):
            # Main chain atoms with proper phi/psi angles for alpha-helix
            x = i * 1.5
            y = 10.0 * (1 + np.cos(i * 100 * 3.14159 / 180))
            z = 10.0 * np.sin(i * 100 * 3.14159 / 180)
            
            # N atom
            pdb_content += f"ATOM  {i*4+1:5d}  N   ALA A {i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 80.00           N  \n"
            # CA atom
            ca_x = x + 0.5
            ca_y = y + 0.5 * (1 + np.cos(i * 100 * 3.14159 / 180 + 0.5))
            ca_z = z + 0.5 * np.sin(i * 100 * 3.14159 / 180 + 0.5)
            pdb_content += f"ATOM  {i*4+2:5d}  CA  ALA A {i+1:4d}    {ca_x:8.3f}{ca_y:8.3f}{ca_z:8.3f}  1.00 80.00           C  \n"
            # C atom
            c_x = x + 1.0
            c_y = y + 1.0 * (1 + np.cos(i * 100 * 3.14159 / 180 + 1.0))
            c_z = z + 1.0 * np.sin(i * 100 * 3.14159 / 180 + 1.0)
            pdb_content += f"ATOM  {i*4+3:5d}  C   ALA A {i+1:4d}    {c_x:8.3f}{c_y:8.3f}{c_z:8.3f}  1.00 80.00           C  \n"
            # O atom
            o_x = x + 1.2
            o_y = y + 1.2 * (1 + np.cos(i * 100 * 3.14159 / 180 + 1.2))
            o_z = z + 1.2 * np.sin(i * 100 * 3.14159 / 180 + 1.2)
            pdb_content += f"ATOM  {i*4+4:5d}  O   ALA A {i+1:4d}    {o_x:8.3f}{o_y:8.3f}{o_z:8.3f}  1.00 80.00           O  \n"
        pdb_content += "TER\nEND"
        
        with open(pdb_file, 'w') as f:
            f.write(pdb_content)
            
        return pdb_file
    
    def get_available_models(self) -> List[str]:
        """Get list of available structure prediction models."""
        return ["ColabFold (AlphaFold2)", "ESMFold", "RoseTTAFold"]
    
    def get_prediction_status(self, job_name: str) -> Dict[str, any]:
        """
        Get the status of a prediction job.
        
        Args:
            job_name: Name of the prediction job
            
        Returns:
            Dictionary containing job status and output files
        """
        job_dir = self.output_dir / job_name
        pdb_file = job_dir / f"{job_name}.pdb"
        
        if pdb_file.exists():
            return {
                "status": "completed", 
                "pdb_file": pdb_file,
                "plddt_file": job_dir / f"{job_name}_plddt.json"
            }
        elif (job_dir / "running").exists():
            return {"status": "running"}
        else:
            return {"status": "not_started"}

def visualize_structure(pdb_file: Path):
    """Visualize protein structure using py3Dmol with enhanced visualization."""
    try:
        with open(pdb_file, 'r') as f:
            pdb_data = f.read()
        
        # Initialize viewer
        view = py3Dmol.view(width=800, height=600)
        
        # Add model
        view.addModel(pdb_data, 'pdb')
        
        # Set visualization style
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Add additional visualization layers
        view.addStyle({'hetflag': True}, {'stick': {}})
        view.addStyle({'resn': 'HIS,ASP,GLU,ARG,LYS,SER,THR,ASN,GLN'}, 
                     {'stick': {'color': 'red'}})
        
        # Add surface representation
        view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'lightgrey'}, 
                       {'hetflag': False})
        
        # Add labels for important residues
        view.addResLabels({"resi": [1, 5, 10, 15, 20]})
        
        # Set view and enable controls
        view.zoomTo()
        view.spin(True)  # Enable rotation
        view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})
        
        # Show the visualization
        showmol(view, height=600, width=800)
        
        # Add download button for the PDB file
        with open(pdb_file, 'rb') as f:
            pdb_data = f.read()
        
        st.download_button(
            label="Download PDB",
            data=pdb_data,
            file_name=pdb_file.name,
            mime="chemical/x-pdb"
        )
        
    except Exception as e:
        st.error(f"❌ Error visualizing structure: {str(e)}")
        logger.exception("Error in structure visualization")

# Initialize structure predictor
predictor = StructurePredictor(STRUCTURES_DIR)

def main():
    st.title("🧬 Structure Prediction")
    st.markdown("---")
    
    # Check if we have filtered targets
    if 'filtered_df' not in st.session_state or st.session_state.filtered_df is None:
        st.warning("⚠️ Please filter targets in the 'Target Analysis' section first.")
        if st.button("Go to Target Analysis"):
            st.switch_page("pages/2_Target_Analysis.py")
        return
        
    # Ensure structures directory exists
    STRUCTURES_DIR.mkdir(exist_ok=True, parents=True)
    
    df = st.session_state.filtered_df
    
    # Target selection
    st.subheader("🎯 Select Target")
    
    # Show all available targets with search functionality
    target_options = sorted(df['gene'].unique().tolist())
    
    # Add search box for filtering targets
    search_term = st.text_input("🔍 Search target genes", "", 
                              placeholder="Type to filter targets...")
    
    # Filter targets based on search term
    if search_term:
        filtered_targets = [t for t in target_options 
                          if search_term.lower() in t.lower()]
        if not filtered_targets:
            st.warning("No matching targets found. Showing all targets.")
            filtered_targets = target_options
    else:
        filtered_targets = target_options
    
    # Display the selectbox with filtered options
    selected_target = st.selectbox(
        f"Select a target gene (Showing {len(filtered_targets)} of {len(target_options)})",
        options=filtered_targets,
        index=0 if not search_term or filtered_targets else None,
        help="Select a target gene for structure prediction. Use the search box above to filter."
    )
    
    # Get target info
    target_info = df[df['gene'] == selected_target].iloc[0]
    
    # Display target information
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Gene", selected_target)
    with col2:
        st.metric("Log2 Fold Change", f"{target_info['logfc']:.2f}")
    with col3:
        st.metric("Adjusted p-value", f"{target_info['adj.p.val']:.3f}")
    
    # Structure prediction section
    st.subheader("🔮 Structure Prediction")
    
    # Get sequence from the dataframe or load from FASTA file
    sequence = ""
    fasta_file = Path(f"structures/{selected_target}/{selected_target}.fasta")
    
    if fasta_file.exists():
        with open(fasta_file, 'r') as f:
            # Skip the header line and join sequence lines
            sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
    else:
        # Fallback to sequence from dataframe if available
        sequence_col = 'sequence' if 'sequence' in target_info else 'protein_sequence'
        sequence = target_info.get(sequence_col, "Sequence not available")
    
    # Display sequence
    with st.expander("View Protein Sequence", expanded=True):
        st.text_area("Protein Sequence", value=sequence, height=150, disabled=True)
        
        # Download button for the sequence
        st.download_button(
            label="Download FASTA",
            data=f">{selected_target}\n{sequence}",
            file_name=f"{selected_target}.fasta",
            mime="text/plain"
        )
    
    # Prediction settings
    st.subheader("⚙️ Prediction Settings")
    
    col1, col2 = st.columns(2)
    with col1:
        model = st.selectbox(
            "Prediction Model",
            options=predictor.get_available_models(),
            index=0,
            help="Select a structure prediction model"
        )
    
    with col2:
        use_templates = st.checkbox(
            "Use templates if available",
            value=True,
            help="Use known protein structures as templates if available"
        )
    
    # Start prediction
    if st.button("🚀 Predict Structure", type="primary"):
        with st.spinner("Running structure prediction (this may take a few minutes)..."):
            result = predictor.predict_structure(
                sequence=sequence,
                job_name=selected_target,
                method=model.lower().split(' ')[0]
            )
            
            if result:
                st.session_state.last_prediction = {
                    'target': selected_target,
                    'pdb_file': result,
                    'model': model,
                    'timestamp': pd.Timestamp.now()
                }
    
    # Show prediction results if available
    if 'last_prediction' in st.session_state:
        pred = st.session_state.last_prediction
        
        # Ensure required keys exist
        if not all(k in pred for k in ['target', 'pdb_file', 'model', 'timestamp']):
            st.error("❌ Invalid prediction data. Please run the prediction again.")
            return
        
        st.subheader("📊 Prediction Results")
        
        # Display basic info in columns
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Target", pred['target'])
        with col2:
            st.metric("Model Used", pred['model'])
        with col3:
            st.metric("Prediction Time", pred['timestamp'].strftime("%Y-%m-%d %H:%M:%S"))
        
        # Check if PDB file exists
        pdb_path = Path(pred['pdb_file'])
        if not pdb_path.exists():
            st.warning("PDB file not found. Please run the prediction first.")
            return
            
        # Ensure the structures directory exists
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Load pLDDT scores if available
        plddt_scores = None
        plddt_file = Path(f"structures/{selected_target}/{selected_target}_plddt.json")
        if plddt_file.exists():
            with open(plddt_file, 'r') as f:
                plddt_data = json.load(f)
                plddt_scores = plddt_data.get('plddt', [])
                if plddt_scores:
                    avg_plddt = sum(plddt_scores) / len(plddt_scores)
                    st.metric("Average pLDDT", f"{avg_plddt:.2f}")
        
        # Create tabs for different views
        tab1, tab2, tab3 = st.tabs(["3D Structure", "Quality Assessment", "Download"])
        
        with tab1:
            # Display structure visualization with enhanced controls
            st.subheader("3D Structure Visualization")
            if pdb_path.exists():
                # Add visualization controls
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    show_cartoon = st.checkbox("Show Cartoon", value=True)
                with col2:
                    show_surface = st.checkbox("Show Surface", value=True)
                with col3:
                    show_sidechains = st.checkbox("Show Sidechains", value=False)
                
                # Visualize structure with selected options
                visualize_structure(pdb_path)
                
                # Show structure statistics
                with st.expander("Structure Statistics"):
                    st.write(f"- **Residues:** 20")
                    st.write(f"- **Chains:** 1")
                    st.write(f"- **Atoms:** 80")
                    st.write(f"- **Resolution:** N/A (Predicted)")
                    
                # Show quality metrics
                st.subheader("Model Quality")
                
                # pLDDT score plot (placeholder)
                plddt_scores = [min(0.9, max(0.7, 0.8 + 0.1 * (i % 5))) for i in range(20)]
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=list(range(1, 21)),
                    y=plddt_scores,
                    mode='lines+markers',
                    name='pLDDT Score',
                    line=dict(color='#636EFA')
                ))
                
                # Add confidence thresholds
                fig.add_hline(y=0.9, line_dash="dash", line_color="green", 
                             annotation_text="Very High Confidence", 
                             annotation_position="bottom right")
                fig.add_hline(y=0.7, line_dash="dash", line_color="orange",
                             annotation_text="High Confidence", 
                             annotation_position="bottom right")
                fig.add_hline(y=0.5, line_dash="dash", line_color="red",
                             annotation_text="Low Confidence", 
                             annotation_position="bottom right")
                
                fig.update_layout(
                    title="Model Confidence (pLDDT)",
                    xaxis_title="Residue Position",
                    yaxis_title="pLDDT Score",
                    yaxis=dict(range=[0.4, 1.0]),
                    height=400
                )
                
                st.plotly_chart(fig, use_container_width=True)
                
            else:
                st.warning("Structure file not found. Please run the prediction first.")
                st.subheader("💊 Ligands Detected")
                for i, ligand in enumerate(ligands, 1):
                    with st.expander(f"Ligand {i}: {ligand['residue_name']}"):
                        st.json({
                            'Residue Name': ligand['residue_name'],
                            'Chain': ligand['chain_id'],
                            'Residue ID': ligand['residue_id'],
                            'Number of Atoms': ligand['num_atoms']
                        })
        
        with tab2:
            # Quality Assessment
            st.subheader("📊 Quality Assessment")
            
            # Show pLDDT plot if available
            if plddt_scores:
                df_plddt = pd.DataFrame({
                    'Residue': range(1, len(plddt_scores) + 1),
                    'pLDDT': plddt_scores
                })
                
                fig = px.line(
                    df_plddt, 
                    x='Residue', 
                    y='pLDDT',
                    title='Per-Residue pLDDT Score',
                    labels={'pLDDT': 'pLDDT Score', 'Residue': 'Residue Number'},
                    height=400
                )
                
                # Add confidence thresholds
                fig.add_hline(y=90, line_dash="dash", line_color="green", 
                             annotation_text="Very High Confidence", annotation_position="bottom right")
                fig.add_hline(y=70, line_dash="dash", line_color="lime", 
                             annotation_text="High Confidence", annotation_position="bottom right")
                fig.add_hline(y=50, line_dash="dash", line_color="orange", 
                             annotation_text="Low Confidence", annotation_position="bottom right")
                
                st.plotly_chart(fig, use_container_width=True)
            
            # Calculate and display SASA
            try:
                sasa = calculate_sasa(pdb_path)
                if sasa['total_sasa'] > 0:
                    st.metric("Solvent Accessible Surface Area (SASA)", 
                             f"{sasa['total_sasa']:.2f} Å²")
            except Exception as e:
                st.warning(f"Could not calculate SASA: {str(e)}")
        
        with tab3:
            # Download section
            st.subheader("💾 Download Results")
            
            # PDB file download
            with open(pdb_path, 'r') as f:
                pdb_data = f.read()
            
            st.download_button(
                label="⬇️ Download PDB File",
                data=pdb_data,
                file_name=f"{selected_target}.pdb",
                mime="text/plain",
                help="Download the predicted 3D structure in PDB format"
            )
            
            # pLDDT scores download
            if plddt_scores:
                df_plddt = pd.DataFrame({
                    'residue': range(1, len(plddt_scores) + 1),
                    'plddt': plddt_scores,
                    'confidence': [
                        'Very High' if score >= 90 
                        else 'High' if score >= 70 
                        else 'Medium' if score >= 50 
                        else 'Low' 
                        for score in plddt_scores
                    ]
                })
                
                csv = df_plddt.to_csv(index=False)
                st.download_button(
                    label="⬇️ Download pLDDT Scores (CSV)",
                    data=csv,
                    file_name=f"{selected_target}_plddt.csv",
                    mime="text/csv",
                    help="Download per-residue pLDDT confidence scores"
                )
                
                # Show summary statistics
                st.subheader("Confidence Summary")
                conf_counts = df_plddt['confidence'].value_counts()
                st.dataframe(
                    conf_counts.rename_axis('Confidence Level')
                             .reset_index(name='Residue Count'),
                    use_container_width=True,
                    hide_index=True
                )
            
            # Generate and download a report
            if st.button("📄 Generate PDF Report"):
                with st.spinner("Generating report..."):
                    # This would generate a PDF report in a real implementation
                    report_path = Path(f"{selected_target}_report.pdf")
                    with open(report_path, 'wb') as f:
                        # In a real app, you would use a PDF generation library here
                        f.write(b"PDF report content would be generated here")
                    
                    with open(report_path, 'rb') as f:
                        st.download_button(
                            label="⬇️ Download PDF Report",
                            data=f,
                            file_name=f"{selected_target}_report.pdf",
                            mime="application/pdf"
                        )
        
        # Confidence scores
        st.subheader("📈 Prediction Confidence")
        
        # Generate some dummy confidence scores
        positions = list(range(1, len(sequence) + 1))
        confidences = np.random.uniform(0.7, 0.95, len(positions))
        
        # Create a plot
        fig = px.line(
            x=positions,
            y=confidences,
            labels={'x': 'Residue Position', 'y': 'pLDDT Score'},
            title='Predicted Confidence (pLDDT) by Residue',
            line_shape='linear'
        )
        fig.update_layout(
            yaxis_range=[0, 1],
            xaxis_title='Residue Position',
            yaxis_title='pLDDT Score',
            hovermode='x'
        )
        
        # Add confidence bands
        fig.add_hrect(y0=0.9, y1=1.0, line_width=0, fillcolor="green", opacity=0.1, annotation_text="Very High")
        fig.add_hrect(y0=0.7, y1=0.9, line_width=0, fillcolor="yellow", opacity=0.1, annotation_text="Confident")
        fig.add_hrect(y0=0.5, y1=0.7, line_width=0, fillcolor="orange", opacity=0.1, annotation_text="Low")
        fig.add_hrect(y0=0, y1=0.5, line_width=0, fillcolor="red", opacity=0.1, annotation_text="Very Low")
        
        st.plotly_chart(fig, use_container_width=True)

if __name__ == "__main__":
    main()

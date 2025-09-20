"""
Nanocarrier Design Page

This module provides functionality for designing and evaluating nanocarriers
for drug delivery based on the properties of the target and ligands.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import yaml
import json
import logging
import base64
from typing import Dict, List, Optional, Tuple
import py3Dmol
from stmol import showmol
from PIL import Image
import io

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set page config
st.set_page_config(
    page_title="Nanocarrier Design - AI Drug Target Platform",
    page_icon="🚀"
)

# Check if we have docking results
if 'docking_results' not in st.session_state:
    st.warning("⚠️ Please run ligand docking first.")
    if st.button("Go to Ligand Docking"):
        st.switch_page("pages/5_Ligand_Docking.py")
    st.stop()

def create_nanocarrier_3d_view(carrier_type: str):
    """Generate a 3D visualization of the nanocarrier."""
    try:
        # This is a simplified representation
        view = py3Dmol.view(width=400, height=300)
        
        if carrier_type == "Liposomes":
            # Simple liposome representation
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 10, 'color': 'lightblue', 'alpha': 0.6})
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 8, 'color': 'white', 'alpha': 0.8})
        elif carrier_type == "Polymeric Nanoparticles":
            # Polymer nanoparticle representation
            for i in range(10):
                x, y, z = np.random.normal(0, 2, 3)
                view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': 1, 'color': 'green', 'alpha': 0.7})
        elif carrier_type == "Dendrimers":
            # Dendrimer representation
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 3, 'color': 'purple', 'alpha': 0.8})
            for i in range(8):
                x, y, z = np.random.normal(0, 1.5, 3)
                view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': 1, 'color': 'orange', 'alpha': 0.7})
        elif carrier_type == "Micelles":
            # Micelle representation
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 8, 'color': 'yellow', 'alpha': 0.6})
            for i in range(12):
                x, y, z = np.random.normal(0, 1, 3)
                view.addSphere({'center': {'x': x*2, 'y': y*2, 'z': z*2}, 'radius': 1, 'color': 'orange', 'alpha': 0.7})
        elif carrier_type == "Gold Nanoparticles":
            # Gold nanoparticle representation
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 5, 'color': 'gold', 'alpha': 0.9})
        elif carrier_type == "Silica Nanoparticles":
            # Silica nanoparticle representation
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 6, 'color': 'lightgray', 'alpha': 0.7})
            view.addSphere({'center': {'x': 0, 'y': 0, 'z': 0}, 'radius': 5, 'color': 'white', 'alpha': 0.8})
        
        view.setStyle({'sphere': {}})
        view.zoomTo()
        view.spin(True)
        return view
    except Exception as e:
        logger.error(f"Error creating 3D view: {str(e)}")
        return None

class NanocarrierDesigner:
    """Class for designing and evaluating nanocarriers for drug delivery."""
    
    def __init__(self):
        self.available_materials = [
            "Liposomes", "Polymeric Nanoparticles", "Dendrimers", 
            "Micelles", "Gold Nanoparticles", "Silica Nanoparticles"
        ]
        
        # Load nanocarrier properties from YAML (in a real app, this would be a file)
        self.carrier_properties = {
            "Liposomes": {
                "size": "50-200 nm",
                "loading_capacity": "High",
                "release_kinetics": "Sustained",
                "stability": "Moderate",
                "surface_modification": "Easy",
                "biocompatibility": "Excellent",
                "cost": "$$",
                "description": "Spherical vesicles with a hydrophilic core and lipid bilayer, ideal for both hydrophobic and hydrophilic drugs."
            },
            "Polymeric Nanoparticles": {
                "size": "20-200 nm",
                "loading_capacity": "High",
                "release_kinetics": "Controlled",
                "stability": "High",
                "surface_modification": "Moderate",
                "biocompatibility": "Good",
                "cost": "$$",
                "description": "Biodegradable polymers that can encapsulate drugs and provide controlled release."
            },
            "Dendrimers": {
                "size": "1-15 nm",
                "loading_capacity": "Moderate",
                "release_kinetics": "Controlled",
                "stability": "High",
                "surface_modification": "Excellent",
                "biocompatibility": "Variable",
                "cost": "$$$",
                "description": "Highly branched, monodisperse macromolecules with multiple functional groups for drug conjugation."
            },
            "Micelles": {
                "size": "10-100 nm",
                "loading_capacity": "Moderate",
                "release_kinetics": "Fast",
                "stability": "Low",
                "surface_modification": "Moderate",
                "biocompatibility": "Good",
                "cost": "$",
                "description": "Amphiphilic molecules that self-assemble in aqueous solutions, suitable for hydrophobic drugs."
            },
            "Gold Nanoparticles": {
                "size": "5-100 nm",
                "loading_capacity": "Low",
                "release_kinetics": "Controlled",
                "stability": "Excellent",
                "surface_modification": "Excellent",
                "biocompatibility": "Good",
                "cost": "$$$",
                "description": "Inert metal nanoparticles with unique optical properties, suitable for theranostic applications."
            },
            "Silica Nanoparticles": {
                "size": "20-200 nm",
                "loading_capacity": "High",
                "release_kinetics": "Controlled",
                "stability": "Excellent",
                "surface_modification": "Excellent",
                "biocompatibility": "Good",
                "cost": "$$",
                "description": "Porous nanoparticles with high surface area for drug loading and controlled release."
            }
        }
        
    def design_nanocarrier(self, carrier_type: str, target_properties: Dict, 
                          ligand_properties: Dict) -> Dict:
        """Design a nanocarrier based on target and ligand properties."""
        # In a real implementation, this would use the properties to design the carrier
        # For now, we'll return some dummy data with enhanced properties
        
        # Base properties based on carrier type
        base_props = {
            "Liposomes": {"size_range": (50, 200), "loading": (60, 90), "stability": (70, 90)},
            "Polymeric Nanoparticles": {"size_range": (20, 200), "loading": (50, 85), "stability": (80, 95)},
            "Dendrimers": {"size_range": (5, 15), "loading": (40, 70), "stability": (90, 100)},
            "Micelles": {"size_range": (10, 100), "loading": (30, 60), "stability": (60, 80)},
            "Gold Nanoparticles": {"size_range": (5, 100), "loading": (20, 50), "stability": (95, 100)},
            "Silica Nanoparticles": {"size_range": (20, 200), "loading": (60, 90), "stability": (85, 100)}
        }
        
        props = base_props.get(carrier_type, {"size_range": (10, 200), "loading": (30, 80), "stability": (70, 100)})
        
        design = {
            'carrier_type': carrier_type,
            'size': np.random.uniform(*props["size_range"]),
            'loading_efficiency': np.random.uniform(*props["loading"]),
            'release_profile': {
                'initial_burst': np.random.uniform(10, 30),
                'sustained_release': np.random.uniform(60, 90, 24).tolist()  # 24h profile
            },
            'stability': np.random.uniform(*props["stability"]),
            'surface_charge': np.random.uniform(-30, 30),
            'polydispersity': np.random.uniform(0.05, 0.3),
            'recommendations': [
                f"Consider {carrier_type}-specific surface modifications",
                "Optimize drug loading parameters for improved efficiency",
                "Test in vitro release profile under physiological conditions"
            ]
        }
        
        # Add targeting ligands if specified
        if ligand_properties and 'name' in ligand_properties:
            design['targeting_ligand'] = ligand_properties['name']
            design['recommendations'].append(
                f"Verify targeting efficiency with {ligand_properties['name']} in vitro"
            )
        
        return design
    
    def recommend_carriers(self, ligand_props: Dict) -> List[Dict]:
        """Recommend nanocarriers based on ligand properties."""
        # In a real implementation, this would use ML or rule-based matching
        # For now, we'll return all carriers with a recommendation score
        
        recommendations = []
        
        for material in self.available_materials:
            # Simple scoring based on ligand properties and carrier characteristics
            score = 0.5  # Base score
            
            # Adjust score based on ligand properties
            if ligand_props['logp'] > 3:  # Hydrophobic ligand
                if material in ["Liposomes", "Micelles"]:
                    score += 0.3
            else:  # Hydrophilic ligand
                if material in ["Dendrimers", "Silica Nanoparticles"]:
                    score += 0.3
            
            # Add some randomness for demo purposes
            score += np.random.uniform(-0.1, 0.1)
            score = max(0.1, min(1.0, score))  # Clamp between 0.1 and 1.0
            
            recommendations.append({
                'material': material,
                'score': score,
                'properties': self.carrier_properties[material]
            })
        
        # Sort by score (descending)
        recommendations.sort(key=lambda x: x['score'], reverse=True)
        
        return recommendations
    
    def get_carrier_comparison_chart(self, recommendations: List[Dict]) -> go.Figure:
        """Generate a radar chart comparing different nanocarriers."""
        categories = ['Loading', 'Stability', 'Release Control', 'Modifiability', 'Biocompatibility']
        
        # Map properties to radar chart values
        property_map = {
            'Loading': 'loading_capacity',
            'Stability': 'stability',
            'Release Control': 'release_kinetics',
            'Modifiability': 'surface_modification',
            'Biocompatibility': 'biocompatibility'
        }
        
        # Convert qualitative properties to numerical scores
        qual_to_num = {
            'Low': 1, 'Moderate': 2, 'High': 3, 'Excellent': 4,
            'Fast': 1, 'Sustained': 2, 'Controlled': 3,
            'Easy': 3, 'Moderate': 2, 'Difficult': 1
        }
        
        # Prepare data for the radar chart
        fig = go.Figure()
        
        # Add a trace for each recommended carrier (top 3)
        for i, rec in enumerate(recommendations[:3]):
            material = rec['material']
            props = rec['properties']
            
            # Convert qualitative properties to numerical scores
            values = [
                qual_to_num.get(props[property_map[cat]], 0) 
                for cat in categories
            ]
            
            fig.add_trace(go.Scatterpolar(
                r=values + [values[0]],  # Close the polygon
                theta=categories + [categories[0]],
                name=material,
                line=dict(color=px.colors.qualitative.Plotly[i]),
                hoverinfo='text',
                hovertext=f"<b>{material}</b><br>" +
                         f"Score: {rec['score']:.2f}<br>" +
                         f"Size: {props['size']}<br>" +
                         f"Cost: {props['cost']}"
            ))
        
        # Update layout
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 4],
                    tickvals=[1, 2, 3, 4],
                    ticktext=['Low', 'Moderate', 'High', 'Excellent']
                )
            ),
            title="Nanocarrier Comparison",
            showlegend=True,
            height=600,
            margin=dict(l=50, r=50, t=80, b=50)
        )
        
        return fig

# Initialize nanocarrier designer
designer = NanocarrierDesigner()

def main():
    st.title("🚀 Nanocarrier Design")
    st.markdown("---")
    
    # Get the docking results
    docking_data = st.session_state.docking_results
    
    # Display target and ligand info
    st.subheader("🎯 Target & Ligand Information")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Target", docking_data['target'])
        st.metric("Docking Program", docking_data['program'])
    with col2:
        # Get the top ligand
        top_ligand = min(docking_data['results'], key=lambda x: x['score'])
        st.metric("Top Ligand", top_ligand['ligand_name'])
        st.metric("Docking Score", f"{top_ligand['score']:.2f} kcal/mol")
    
    # Ligand properties (in a real app, these would be calculated)
    st.subheader("🧪 Ligand Properties")
    
    # Generate some dummy ligand properties
    ligand_props = {
        'mw': np.random.uniform(200, 600),
        'logp': np.random.uniform(-2, 5),
        'hbd': np.random.randint(0, 5),
        'hba': np.random.randint(1, 8),
        'tpsa': np.random.uniform(20, 120),
        'rotatable_bonds': np.random.randint(0, 10),
        'aromatic_rings': np.random.randint(1, 4)
    }
    
    # Display ligand properties
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Molecular Weight", f"{ligand_props['mw']:.1f} g/mol")
        st.metric("H-Bond Donors", ligand_props['hbd'])
    with col2:
        st.metric("LogP", f"{ligand_props['logp']:.2f}")
        st.metric("H-Bond Acceptors", ligand_props['hba'])
    with col3:
        st.metric("TPSA", f"{ligand_props['tpsa']:.1f} Å²")
        st.metric("Rotatable Bonds", ligand_props['rotatable_bonds'])
    with col4:
        st.metric("Aromatic Rings", ligand_props['aromatic_rings'])
    
    # Nanocarrier design
    st.subheader("🧬 Nanocarrier Recommendations")
    
    # Get recommendations
    recommendations = designer.recommend_carriers(ligand_props)
    
    # Display top recommendations
    st.markdown("### Top Recommended Nanocarriers")
    
    # Show top 3 recommendations as cards
    cols = st.columns(3)
    for i, rec in enumerate(recommendations[:3]):
        with cols[i]:
            with st.expander(f"**{i+1}. {rec['material']}** (Score: {rec['score']:.2f})", expanded=True):
                st.markdown(f"**Description:** {rec['properties']['description']}")
                st.markdown(f"**Size:** {rec['properties']['size']}")
                st.markdown(f"**Loading Capacity:** {rec['properties']['loading_capacity']}")
                st.markdown(f"**Release Kinetics:** {rec['properties']['release_kinetics']}")
                st.markdown(f"**Stability:** {rec['properties']['stability']}")
                st.markdown(f"**Surface Modification:** {rec['properties']['surface_modification']}")
                st.markdown(f"**Biocompatibility:** {rec['properties']['biocompatibility']}")
                st.markdown(f"**Cost:** {rec['properties']['cost']}")
    
    # Show comparison chart
    st.subheader("📊 Nanocarrier Comparison")
    
    fig = designer.get_carrier_comparison_chart(recommendations)
    st.plotly_chart(fig, use_container_width=True)
    
    # Custom nanocarrier design
    st.subheader("🎨 Custom Nanocarrier Design")
    
    # Initialize session state for nanocarrier design if it doesn't exist
    if 'nanocarrier_design' not in st.session_state:
        st.session_state.nanocarrier_design = None
    
    with st.form("custom_carrier_form"):
        st.markdown("### Design Your Own Nanocarrier")
        
        col1, col2 = st.columns(2)
        
        with col1:
            material = st.selectbox(
                "Base Material",
                options=designer.available_materials,
                index=0
            )
            
            size = st.slider(
                "Target Size (nm)",
                min_value=10,
                max_value=500,
                value=100,
                step=10
            )
            
            surface_charge = st.select_slider(
                "Surface Charge",
                options=['Strongly Negative', 'Weakly Negative', 'Neutral', 'Weakly Positive', 'Strongly Positive'],
                value='Neutral'
            )
        
        with col2:
            targeting_ligand = st.multiselect(
                "Targeting Ligands",
                options=["None", "Antibodies", "Aptamers", "Peptides", "Folic Acid", "Transferrin"],
                default=["None"]
            )
            
            release_trigger = st.selectbox(
                "Release Trigger",
                options=["pH", "Enzyme", "Redox", "Temperature", "Light", "Ultrasound"],
                index=0
            )
            
            coating = st.selectbox(
                "Surface Coating",
                options=["None", "PEG", "Chitosan", "Hyaluronic Acid", "Polysorbate 80"],
                index=0
            )
        
        # Submit button
        submitted = st.form_submit_button("Design Nanocarrier")
        
        if submitted:
            st.session_state.nanocarrier_design = {
                'material': material,
                'size': size,
                'surface_charge': surface_charge,
                'targeting_ligands': targeting_ligand,
                'release_trigger': release_trigger,
                'coating': coating,
                'ligand': top_ligand['ligand_name'],
                'target': docking_data['target']
            }
    
    # Display the design and download button outside the form
    if st.session_state.nanocarrier_design:
        design = st.session_state.nanocarrier_design
        
        st.success("🎉 Nanocarrier design ready!")
        
        # Show a summary of the design
        with st.expander("📝 Design Summary", expanded=True):
            st.markdown("### Custom Nanocarrier Design")
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown(f"**Base Material:** {design['material']}")
                st.markdown(f"**Size:** {design['size']} nm")
                st.markdown(f"**Surface Charge:** {design['surface_charge']}")
            
            with col2:
                st.markdown(f"**Targeting Ligands:** {', '.join(design['targeting_ligands']) if design['targeting_ligands'] else 'None'}")
                st.markdown(f"**Release Trigger:** {design['release_trigger']}")
                st.markdown(f"**Surface Coating:** {design['coating']}")
        
        st.markdown("### 🎨 Nanocarrier Visualization")
        
        # Local fallback images (base64 encoded to avoid external dependencies)
        local_images = {
            "Liposomes": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iODAiIGZpbGw9IiNhY2U0ZTkiIHN0cm9rZT0iIzY5YjJkZSIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMjAwIiBjeT0iMTUwIiByPSI2MCIgZmlsbD0id2hpdGUiLz48dGV4dCB4PSIyMDAiIHk9IjE1MCIgZm9udC1mYW1pbHk9IkFyaWFsIiBmb250LXNpemU9IjEyIiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBhbGlnbm1lbnQtYmFzZWxpbmU9Im1pZGRsZSIgZmlsbD0iIzMzMyI+TGlwb3NvbWU8L3RleHQ+PC9zdmc+",
            "Polymeric Nanoparticles": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iODAiIGZpbGw9IiNhZGQ4YzciIHN0cm9rZT0iIzhjYjNhNiIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMTUwIiBjeT0iMTUwIiByPSIyMCIgZmlsbD0iI2ZmYzEwNyIvPjxjaXJjbGUgY3g9IjI1MCIgY3k9IjE1MCIgcj0iMjAiIGZpbGw9IiNmZmMxMDciLz48Y2lyY2xlIGN4PSIyMDAiIGN5PSIyMDAiIHI9IjIwIiBmaWxsPSIjZmZjMTA3Ii8+PGNpcmNsZSBjeD0iMTUwIiBjeT0iMTAwIiByPSIyMCIgZmlsbD0iI2ZmYzEwNyIvPjx0ZXh0IHg9IjIwMCIgeT0iMTUwIiBmb250LWZhbWlseT0iQXJpYWwiIGZvbnQtc2l6ZT0iMTIiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGFsaWdubWVudC1iYXNlbGluZT0ibWlkZGxlIiBmaWxsPSIjMzMzIj5Qb2x5bWVyaWMgTmFub3BhcnRpY2xlczwvdGV4dD48L3N2Zz4=",
            "Dendrimers": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iNDAiIGZpbGw9IiNlYjBkZmYiIHN0cm9rZT0iI2I4NGQ5YSIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMjAwIiBjeT0iODAiIHI9IjIwIiBmaWxsPSIjZjQ4MWQ2Ii8+PGNpcmNsZSBjeD0iMjAwIiBjeT0iMjIwIiByPSIyMCIgZmlsbD0iI2Y0ODFkNiIvPjxjaXJjbGUgY3g9IjgwIiBjeT0iMTUwIiByPSIyMCIgZmlsbD0iI2Y0ODFkNiIvPjxjaXJjbGUgY3g9IjMyMCIgY3k9IjE1MCIgcj0iMjAiIGZpbGw9IiNmNDgxZDYiLz48dGV4dCB4PSIyMDAiIHk9IjE1MCIgZm9udC1mYW1pbHk9IkFyaWFsIiBmb250LXNpemU9IjEyIiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBhbGlnbm1lbnQtYmFzZWxpbmU9Im1pZGRsZSIgZmlsbD0iIzMzMyI+RGVuZHJpbWVyczwvdGV4dD48L3N2Zz4=",
            "Micelles": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iODAiIGZpbGw9IiNmZmUwYjUiIHN0cm9rZT0iI2ZmYzEwNyIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMTUwIiBjeT0iMTUwIiByPSI2MCIgZmlsbD0id2hpdGUiLz48dGV4dCB4PSIyMDAiIHk9IjE1MCIgZm9udC1mYW1pbHk9IkFyaWFsIiBmb250LXNpemU9IjEyIiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBhbGlnbm1lbnQtYmFzZWxpbmU9Im1pZGRsZSIgZmlsbD0iIzMzMyI+TWljZWxsZXM8L3RleHQ+PC9zdmc+",
            "Gold Nanoparticles": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iODAiIGZpbGw9IiNGRkQ3MDAiIHN0cm9rZT0iI0ZGQjQwMCIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMTUwIiBjeT0iMTUwIiByPSI2MCIgZmlsbD0iI0ZGRUIzQyIvPjx0ZXh0IHg9IjIwMCIgeT0iMTUwIiBmb250LWZhbWlseT0iQXJpYWwiIGZvbnQtc2l6ZT0iMTIiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGFsaWdubWVudC1iYXNlbGluZT0ibWlkZGxlIiBmaWxsPSIjMzMzIj5Hb2xkIE5hbm9wYXJ0aWNsZXM8L3RleHQ+PC9zdmc+",
            "Silica Nanoparticles": "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI0MDAiIGhlaWdodD0iMzAwIiB2aWV3Qm94PSIwIDAgNDAwIDMwMCI+PHJlY3Qgd2lkdGg9IjEwMCUiIGhlaWdodD0iMTAwJSIgZmlsbD0iI2VlZSIvPjxjaXJjbGUgY3g9IjIwMCIgY3k9IjE1MCIgcj0iODAiIGZpbGw9IiNkZGQiIHN0cm9rZT0iI2JiYiIgc3Ryb2tlLXdpZHRoPSIyIi8+PGNpcmNsZSBjeD0iMjAwIiBjeT0iMTUwIiByPSI2MCIgZmlsbD0iI2VlZSIvPjx0ZXh0IHg9IjIwMCIgeT0iMTUwIiBmb250LWZhbWlseT0iQXJpYWwiIGZvbnQtc2l6ZT0iMTIiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGFsaWdubWVudC1iYXNlbGluZT0ibWlkZGxlIiBmaWxsPSIjMzMzIj5TaWxpY2EgTmFub3BhcnRpY2xlczwvdGV4dD48L3N2Zz4="
        }
        
        # Try to show 3D visualization, fall back to static image
        with st.spinner("Generating 3D visualization..."):
            try:
                view = create_nanocarrier_3d_view(design['material'])
                if view:
                    showmol(view, height=400, width=600)
                    st.caption(f"3D visualization of {design['material']}")
                else:
                    raise Exception("3D view creation failed")
            except Exception as e:
                logger.warning(f"3D visualization failed: {str(e)}")
                st.warning("3D visualization not available. Showing static representation instead.")
                # Use the appropriate image based on carrier type
                image_data = local_images.get(
                    design['material'],
                    local_images["Liposomes"]  # Default fallback
                )
                st.image(
                    image_data,
                    use_container_width=True,
                    caption=f"{design['material']} nanocarrier with {', '.join(design['targeting_ligands']) if design['targeting_ligands'] else 'no'} targeting ligands"
                )
        
        # Save design button (outside the form)
        st.download_button(
            label="💾 Save Nanocarrier Design",
            data=json.dumps(design, indent=2),
            file_name=f"nanocarrier_{design['target']}_{design['material'].lower().replace(' ', '_')}.json",
            mime="application/json"
        )

if __name__ == "__main__":
    main()

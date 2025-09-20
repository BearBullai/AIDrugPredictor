"""
Ligand Docking Page

This module provides functionality for docking small molecules into protein binding sites
using various docking programs and scoring functions.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import base64
import logging
import tempfile
from typing import Dict, List, Optional, Tuple
from PIL import Image
import io

# Try optional RDKit imports (not always available in cloud builds)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config is managed in app.py

# Check if we have detected pockets
if 'detected_pockets' not in st.session_state:
    st.warning("⚠️ Please detect binding pockets first.")
    if st.button("Go to Pocket Detection"):
        st.switch_page("pages/4_Pocket_Detection.py")
    st.stop()

class DockingEngine:
    """Class for handling molecular docking operations."""
    
    def __init__(self):
        self.available_programs = ["AutoDock Vina", "GNINA", "QuickVina", "SMINA"]
        self.available_databases = [
            "DrugBank", "ChEMBL", "ZINC15", "Custom"
        ]
        
    def prepare_ligands(self, input_file: Path, format: str = "sdf") -> List[Dict]:
        """Prepare ligands for docking."""
        # In a real implementation, this would convert and prepare ligand files
        # For now, we'll return some dummy data
        
        # Generate some sample ligands
        ligands = []
        num_ligands = 20
        
        for i in range(1, num_ligands + 1):
            ligand = {
                'id': f"LIG{i:03d}",
                'name': f"Compound_{i}",
                'smiles': f"CCOC(=O)C{i}N",
                'mw': np.random.uniform(200, 600),
                'logp': np.random.uniform(-2, 5),
                'hbd': np.random.randint(0, 5),
                'hba': np.random.randint(1, 8),
                'tpsa': np.random.uniform(20, 120)
            }
            ligands.append(ligand)
        
        return ligands
    
    def run_docking(self, protein_file: Path, pocket: Dict, ligands: List[Dict], 
                   program: str = "AutoDock Vina", num_poses: int = 3) -> List[Dict]:
        """Run molecular docking with enhanced results including 2D interaction diagrams."""
        # In a real implementation, this would call the actual docking software
        # For now, we'll return some dummy results with enhanced interaction data
        
        results = []
        
        for i, ligand in enumerate(ligands[:10]):  # Limit to top 10 for demo
            # Generate a simple 2D structure for the ligand if RDKit is available
            mol = None
            if RDKit_AVAILABLE:
                try:
                    mol = Chem.MolFromSmiles(ligand['smiles'])
                    if mol:
                        AllChem.Compute2DCoords(mol)
                        ligand['mol'] = mol
                except Exception:
                    mol = None
            
            for pose in range(1, num_poses + 1):
                # Generate interaction data
                interactions = {
                    'h_bonds': np.random.randint(1, 5),
                    'hydrophobic': np.random.randint(2, 8),
                    'electrostatic': np.random.randint(0, 3),
                    'residues': [
                        f'ASP{np.random.randint(50, 100)}',
                        f'GLU{np.random.randint(50, 100)}',
                        f'LYS{np.random.randint(50, 100)}',
                        f'PHE{np.random.randint(50, 100)}'
                    ][:np.random.randint(2, 5)]
                }
                
                result = {
                    'ligand_id': ligand['id'],
                    'ligand_name': ligand['name'],
                    'smiles': ligand.get('smiles', ''),
                    'mol': ligand.get('mol'),
                    'pose': pose,
                    'score': np.random.uniform(-12.5, -5.0),
                    'rmsd_lb': np.random.uniform(0.0, 2.0),
                    'rmsd_ub': np.random.uniform(0.0, 2.0) + 0.5,
                    'interactions': interactions,
                    'file': f"docking_pose_{ligand['id']}_pose{pose}.pdbqt"
                }
                results.append(result)
        
        # Sort by score (ascending for Vina, lower is better)
        results.sort(key=lambda x: x['score'])
        
        # Add rank
        for i, result in enumerate(results, 1):
            result['rank'] = i
        
        return results
    
    def get_docking_visualization(self, protein_file: Path, results: List[Dict]) -> go.Figure:
        """Generate a 3D visualization of the docking results."""
        # This is a placeholder for actual 3D visualization
        # In a real implementation, you would use NGLView or Py3Dmol
        
        # Create a simple 3D scatter plot with plotly
        fig = go.Figure()
        
        # Add protein backbone (dummy data)
        x = np.linspace(-10, 10, 100)
        y = np.sin(x) * 2
        z = np.cos(x) * 2
        
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=4, color='lightblue'),
            name='Protein Backbone',
            opacity=0.3
        ))
        
        # Add top 5 poses
        top_results = sorted(results, key=lambda x: x['score'])[:5]
        colors = px.colors.qualitative.Plotly
        
        for i, result in enumerate(top_results):
            color = colors[i % len(colors)]
            
            # Generate some random coordinates around the pocket
            x_pos = np.random.normal(0, 3, 20) + np.random.uniform(-2, 2)
            y_pos = np.random.normal(0, 3, 20) + np.random.uniform(-2, 2)
            z_pos = np.random.normal(0, 3, 20) + np.random.uniform(-2, 2)
            
            fig.add_trace(go.Scatter3d(
                x=x_pos,
                y=y_pos,
                z=z_pos,
                mode='markers',
                marker=dict(
                    size=4,
                    color=color,
                    opacity=0.8,
                    line=dict(width=1, color='black')
                ),
                name=f"{result['ligand_name']} (Score: {result['score']:.2f})",
                text=f"Ligand: {result['ligand_name']}\n"
                     f"Pose: {result['pose']}\n"
                     f"Score: {result['score']:.2f} kcal/mol\n"
                     f"H-bonds: {result['interactions']['h_bonds']}",
                hoverinfo='text+name'
            ))
        
        # Update layout
        fig.update_layout(
            title="Top 5 Docking Poses",
            scene=dict(
                xaxis_title='X (Å)',
                yaxis_title='Y (Å)',
                zaxis_title='Z (Å)',
                aspectmode='data'
            ),
            margin=dict(l=0, r=0, b=0, t=30),
            height=600,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=0.01
            )
        )
        
        return fig

# Initialize docking engine
docking_engine = DockingEngine()

def main():
    st.title("⚗️ Ligand Docking")
    st.markdown("---")
    
    # Get the last prediction and pocket detection results
    pred = st.session_state.last_prediction
    pockets_data = st.session_state.detected_pockets
    
    # Display target info
    st.subheader("🎯 Target Information")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Target", pred['target'])
        st.metric("PDB File", Path(pred['pdb_file']).name)
    with col2:
        st.metric("Prediction Model", pred['model'])
        st.metric("Pocket Detection Method", pockets_data['method'])
    
    # Pocket selection
    st.subheader("🕳️ Select Binding Pocket")
    
    # Get pocket options
    pocket_options = [f"Pocket {p['id']} (Score: {p['score']:.2f})" for p in pockets_data['pockets'][:5]]
    selected_pocket = st.selectbox(
        "Select a pocket for docking",
        options=pocket_options,
        index=0,
        help="Select a binding pocket to dock ligands into"
    )
    
    # Get the selected pocket
    pocket_id = int(selected_pocket.split()[1])
    pocket = next(p for p in pockets_data['pockets'] if p['id'] == pocket_id)
    
    # Docking settings
    st.subheader("⚙️ Docking Settings")
    
    col1, col2 = st.columns(2)
    with col1:
        program = st.selectbox(
            "Docking Program",
            options=docking_engine.available_programs,
            index=0,
            help="Select a docking program to use"
        )
        
        database = st.selectbox(
            "Ligand Database",
            options=docking_engine.available_databases,
            index=0,
            help="Select a database of compounds to dock"
        )
    
    with col2:
        num_poses = st.slider(
            "Number of poses per ligand",
            min_value=1,
            max_value=10,
            value=3,
            help="Number of binding poses to generate per ligand"
        )
        
        exhaustiveness = st.slider(
            "Exhaustiveness",
            min_value=1,
            max_value=100,
            value=8,
            help="Exhaustiveness of the search (higher values give better results but take longer)"
        )
    
    # Upload custom ligands if needed
    custom_ligands = None
    if database == "Custom":
        uploaded_files = st.file_uploader(
            "Upload ligand files (SDF, MOL2, PDBQT)",
            type=["sdf", "mol2", "pdbqt", "pdb"],
            accept_multiple_files=True
        )
        
        if uploaded_files:
            custom_ligands = []
            for file in uploaded_files:
                custom_ligands.append({
                    'name': file.name,
                    'content': file.getvalue()
                })
    
    # Run docking
    if st.button("⚡ Run Docking", type="primary"):
        with st.spinner("Running docking simulation..."):
            # Prepare ligands
            if database == "Custom" and custom_ligands:
                # Use custom ligands
                ligands = [{'id': f"CUST_{i}", 'name': lig['name']} for i, lig in enumerate(custom_ligands, 1)]
            else:
                # Use built-in database (in a real app, this would load from the actual database)
                ligands = docking_engine.prepare_ligands(Path("dummy.sdf"))
            
            # Run docking
            results = docking_engine.run_docking(
                protein_file=Path(pred['pdb_file']),
                pocket=pocket,
                ligands=ligands,
                program=program,
                num_poses=num_poses
            )
            
            st.session_state.docking_results = {
                'target': pred['target'],
                'pocket_id': pocket_id,
                'program': program,
                'ligands': ligands,
                'results': results,
                'timestamp': pd.Timestamp.now()
            }
    
    # Show results if available
    if 'docking_results' in st.session_state:
        docking_data = st.session_state.docking_results
        
        st.subheader("📊 Docking Results")
        
        # Show top poses
        st.subheader("🏆 Top Scoring Poses")
        
        # Create a dataframe of results
        results_df = pd.DataFrame([{
            'Rank': r['rank'],
            'Ligand ID': r['ligand_id'],
            'Ligand Name': r['ligand_name'],
            'Pose': r['pose'],
            'Docking Score': f"{r['score']:.2f} kcal/mol",
            'H-Bonds': r['interactions']['h_bonds'],
            'Hydrophobic': r['interactions']['hydrophobic'],
            'Electrostatic': r['interactions']['electrostatic'],
            'RMSD LB': f"{r['rmsd_lb']:.2f}",
            'RMSD UB': f"{r['rmsd_ub']:.2f}"
        } for r in docking_data['results'][:20]])  # Show top 20 results
        
        # Display the results table with proper formatting
        # First, create a copy for display with formatted strings
        display_df = results_df.copy()
        
        # For sorting and styling, we need numeric values
        numeric_df = results_df.copy()
        numeric_df['Docking Score'] = numeric_df['Docking Score'].str.replace(' kcal/mol', '').astype(float)
        
        # Apply styling to the display dataframe
        def color_score(val):
            # Color based on the numeric value (lower is better for docking scores)
            try:
                score = float(str(val).replace(' kcal/mol', ''))
                # Scale from red (bad) to green (good) for scores between -5 and -15
                # More negative (better) scores get greener
                if score < -12:
                    return 'background-color: rgba(0, 200, 0, 0.3)'  # Light green
                elif score < -10:
                    return 'background-color: rgba(200, 200, 0, 0.3)'  # Light yellow
                else:
                    return 'background-color: rgba(200, 0, 0, 0.1)'  # Light red
            except:
                return ''
        
        # Apply styling
        styled_df = display_df.style.applymap(color_score, subset=['Docking Score'])
        
        # Display the styled dataframe
        st.dataframe(
            styled_df,
            height=500,
            use_container_width=True
        )
        
        # 3D Visualization
        st.subheader("🔄 3D Visualization")
        
        # Show the 3D plot
        fig = docking_engine.get_docking_visualization(
            Path(pred['pdb_file']),
            docking_data['results']
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Interaction analysis
        st.subheader("🔍 Interaction Analysis")
        
        # Get the top result
        top_result = min(docking_data['results'], key=lambda x: x['score'])
        
        # Display interaction details
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Top Score", f"{top_result['score']:.2f} kcal/mol")
        with col2:
            st.metric("H-Bonds", top_result['interactions']['h_bonds'])
        with col3:
            st.metric("Hydrophobic Interactions", top_result['interactions']['hydrophobic'])
        
        # Generate 2D interaction diagram (only if RDKit is available)
        try:
            if not RDKit_AVAILABLE:
                st.info("RDKit is not available; rendering 2D diagram with SmilesDrawer.")
                smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Example: Aspirin
                smiles_html = f"""
                <div id=\"mol-canvas\" style=\"width: 100%; text-align: center;\"></div>
                <script src=\"https://unpkg.com/smiles-drawer@2.0.1/dist/smiles-drawer.min.js\"></script>
                <script>
                  const smiles = '{smiles}';
                  const options = {{ width: 400, height: 300, bondThickness: 1.0 }};
                  const viewer = new SmilesDrawer.Drawer(options);
                  SmilesDrawer.parse(smiles, function(tree) {{
                    const svg = viewer.draw(tree, 'light', false);
                    document.getElementById('mol-canvas').innerHTML = svg;
                  }}, function(err) {{
                    document.getElementById('mol-canvas').innerHTML = '<div style="padding: 1rem; color: #b00;">Failed to render molecule.</div>';
                  }});
                </script>
                """
                st.components.v1.html(smiles_html, height=330)
                # Also show interaction details
                st.markdown("### Interaction Details")
                st.markdown(f"- **H-Bonds:** {top_result['interactions']['h_bonds']}")
                st.markdown(f"- **Hydrophobic Interactions:** {top_result['interactions']['hydrophobic']}")
                st.markdown(f"- **Electrostatic Interactions:** {top_result['interactions']['electrostatic']}")
                # Skip RDKit branch
                raise SystemExit
            # Generate a sample molecule (in a real app, this would be the actual ligand)
            smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Example: Aspirin
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is not None:
                # Generate 2D coordinates
                rdDepictor.Compute2DCoords(mol)
                
                # Create drawer
                drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                
                # Convert to image
                svg = drawer.GetDrawingText()
                svg = svg.replace('svg:', '')
                
                # Display the molecule
                st.components.v1.html(f"""
                <div style="width:100%; text-align:center;">
                    {svg}
                    <p>2D structure of {top_result['ligand_name']} (Pose {top_result['pose']})</p>
                </div>
                """, height=350)
                
                # Add interaction details
                st.markdown("### Interaction Details")
                st.markdown(f"- **H-Bonds:** {top_result['interactions']['h_bonds']}")
                st.markdown(f"- **Hydrophobic Interactions:** {top_result['interactions']['hydrophobic']}")
                st.markdown(f"- **Electrostatic Interactions:** {top_result['interactions']['electrostatic']}")
                
            else:
                st.warning("Could not generate 2D structure for the ligand.")
            
        except SystemExit:
            pass
        except Exception as e:
            if RDKit_AVAILABLE:
                st.error(f"Error generating interaction diagram: {str(e)}")
            # Fallback placeholder
            st.markdown("_2D interaction diagram not available._")
        
        # Save results
        st.download_button(
            label="💾 Save Docking Results",
            data=results_df.to_csv(index=False).encode('utf-8'),
            file_name=f"{docking_data['target']}_pocket{docking_data['pocket_id']}_docking.csv",
            mime="text/csv"
        )

if __name__ == "__main__":
    main()

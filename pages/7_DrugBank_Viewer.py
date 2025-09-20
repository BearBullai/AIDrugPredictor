# Try optional RDKit imports (not always available in cloud builds)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False

"""
DrugBank Viewer Page

This module provides functionality for viewing and analyzing DrugBank entries.
"""
import streamlit as st
import pandas as pd
import requests
from pathlib import Path
import json
from typing import Dict, Optional, List
from PIL import Image
import io
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Page config is managed in app.py

class DrugBankViewer:
    """Class for handling DrugBank data retrieval and visualization."""
    
    def __init__(self):
        self.cache_dir = Path("data/drugbank")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def get_drug_info(self, drug_id: str) -> Optional[Dict]:
        """Get drug information from DrugBank."""
        cache_file = self.cache_dir / f"{drug_id}.json"
        
        # Try to load from cache first
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        # In a real implementation, this would make an API call to DrugBank
        # For now, we'll use mock data
        mock_data = self._generate_mock_drug_data(drug_id)
        
        # Save to cache
        with open(cache_file, 'w') as f:
            json.dump(mock_data, f)
        
        return mock_data
    
    def _generate_mock_drug_data(self, drug_id: str) -> Dict:
        """Generate mock drug data for demonstration."""
        # Create a simple molecule for visualization (keep SMILES regardless of RDKit availability)
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin-like structure
        
        return {
            "id": drug_id,
            "name": f"Drug {drug_id}",
            "description": f"Description for drug {drug_id}. This is a placeholder description that would contain detailed information about the drug's properties, mechanism of action, and clinical uses.",
            "indications": ["Pain", "Fever", "Inflammation"],
            "mechanism_of_action": f"Drug {drug_id} works by inhibiting the production of certain natural substances that cause pain, fever, and inflammation.",
            "molecular_formula": "C9H8O4",
            "molecular_weight": 180.16,
            "targets": ["COX-1", "COX-2"],
            "smiles": smiles,
            "status": "Approved",
            "atc_codes": ["N02BA01", "B01AC06"],
            "half_life": "2-3 hours",
            "protein_binding": "80-90%",
            "metabolism": "Hepatic"
        }
    
    def draw_molecule(self, smiles: str, size: tuple = (400, 300)) -> Optional[Image.Image]:
        """Generate a 2D structure image from SMILES with improved error handling."""
        try:
            if not RDKIT_AVAILABLE:
                logger.info("RDKit not available; cannot draw molecule.")
                return None
            # Convert SMILES to molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Could not parse SMILES: {smiles}")
                return None
            
            # Generate 2D coordinates
            mol = Chem.Mol(mol.ToBinary())
            AllChem.Compute2DCoords(mol)
            
            # Try different drawing methods for compatibility
            try:
                # Try using the modern drawing method first
                from rdkit.Chem.Draw import rdMolDraw2D
                
                # Create drawer
                drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
                
                # Set drawing options
                drawer.drawOptions().useBWAtomPalette()
                drawer.drawOptions().addAtomIndices = False
                drawer.drawOptions().addStereoAnnotation = True
                
                # Draw and get the image
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                img_data = drawer.GetDrawingText()
                return Image.open(io.BytesIO(img_data))
                
            except Exception as e:
                logger.warning(f"Modern drawing method failed, falling back to basic method: {str(e)}")
                # Fallback to basic drawing method
                try:
                    from rdkit.Chem import Draw
                    return Draw.MolToImage(mol, size=size)
                except Exception as e2:
                    logger.error(f"Basic drawing method also failed: {str(e2)}")
                    return None
                    
        except Exception as e:
            logger.error(f"Error in draw_molecule: {str(e)}")
            return None
    
    def _show_fallback_image(self, drug_name: str):
        """Display a fallback image when molecule generation fails."""
        fallback_svg = """
        <svg width="100%" height="200" viewBox="0 0 400 200" xmlns="http://www.w3.org/2000/svg">
            <rect width="100%" height="100%" fill="#f5f5f5"/>
            <circle cx="200" cy="100" r="50" fill="#e0e0e0" stroke="#9e9e9e" stroke-width="2"/>
            <text x="200" y="100" font-family="Arial" font-size="14" text-anchor="middle" dominant-baseline="middle" fill="#616161">
                Structure not available
            </text>
            <text x="200" y="130" font-family="Arial" font-size="12" text-anchor="middle" fill="#9e9e9e">
                {drug_name}
            </text>
        </svg>
        """.format(drug_name=drug_name)
        st.components.v1.html(fallback_svg, height=220)

def display_drug_info(drug_id: str):
    """Display drug information in the UI."""
    viewer = DrugBankViewer()
    
    with st.spinner(f"Loading information for {drug_id}..."):
        drug_info = viewer.get_drug_info(drug_id)
    
    if not drug_info:
        st.error(f"❌ Could not find drug with ID: {drug_id}")
        st.info("Please check the DrugBank ID and try again.")
        return
    
    # Display basic info in two columns
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.title(f"💊 {drug_info.get('name', 'Unknown Drug')}")
        
        # Display molecule image with better error handling and loading state
        with st.spinner("Generating molecule visualization..."):
            try:
                if 'smiles' in drug_info and drug_info['smiles']:
                    # Try to generate 2D structure
                    img = viewer.draw_molecule(drug_info['smiles'], size=(400, 300))
                    if img:
                        # Convert to bytes for better display
                        buffered = io.BytesIO()
                        img.save(buffered, format="PNG")
                        st.image(
                            buffered, 
                            use_container_width=True,
                            caption=f"2D Structure of {drug_info.get('name', 'the drug')}"
                        )
                        
                        # Add download button for the structure
                        st.download_button(
                            label="⬇️ Download 2D Structure",
                            data=buffered.getvalue(),
                            file_name=f"{drug_id}_2d_structure.png",
                            mime="image/png"
                        )
                    else:
                        # RDKit path failed; render with SmilesDrawer in-browser
                        if not RDKIT_AVAILABLE:
                            st.info("RDKit is not available; rendering 2D diagram with SmilesDrawer.")
                        smiles = drug_info['smiles']
                        smiles_html = f"""
                        <div style=\"width:100%; text-align:center;\">
                          <canvas id=\"db-mol-canvas\" width=\"400\" height=\"300\"></canvas>
                        </div>
                        <script src=\"https://cdn.jsdelivr.net/npm/smiles-drawer@2.0.1/dist/smiles-drawer.min.js\"></script>
                        <script>
                          (function() {{
                            const smiles = '{smiles}';
                            const options = {{ width: 400, height: 300, bondThickness: 1.0 }};
                            const drawer = new SmilesDrawer.Drawer(options);
                            SmilesDrawer.parse(smiles)
                              .then(tree => drawer.draw(tree, 'db-mol-canvas', 'light', false))
                              .catch(err => {{
                                const ctx = document.getElementById('db-mol-canvas').getContext('2d');
                                ctx.fillStyle = '#b00';
                                ctx.font = '14px Arial';
                                ctx.fillText('Failed to render molecule', 60, 150);
                              }});
                          }})();
                        </script>
                        """
                        st.components.v1.html(smiles_html, height=330)
                else:
                    st.warning("No SMILES structure available for this drug.")
                    viewer._show_fallback_image(drug_info.get('name', 'the drug'))
                        
            except Exception as e:
                logger.error(f"Error generating molecule image: {str(e)}")
                st.error("Failed to generate molecule visualization. Please try again later.")
                viewer._show_fallback_image(drug_info.get('name', 'the drug'))
    
        # Basic properties
        st.metric("DrugBank ID", drug_id)
        st.metric("Status", drug_info.get("status", "Unknown"))
        st.metric("Molecular Weight", f"{drug_info.get('molecular_weight', 0):.2f} g/mol")
        st.metric("Molecular Formula", drug_info.get("molecular_formula", "N/A"))
    
    with col2:
        st.markdown("### Description")
        st.write(drug_info.get("description", "No description available."))
        
        st.markdown("### Indications")
        for indication in drug_info.get("indications", ["No indications available."]):
            st.markdown(f"- {indication}")
        
        st.markdown("### Mechanism of Action")
        st.write(drug_info.get("mechanism_of_action", "Not available."))
        
        st.markdown("### Targets")
        for target in drug_info.get("targets", ["No targets available."]):
            st.markdown(f"- {target}")
        
        # Additional metadata
        with st.expander("Additional Information"):
            st.markdown("#### Pharmacokinetics")
            st.markdown(f"- **Half-life:** {drug_info.get('half_life', 'N/A')}")
            st.markdown(f"- **Protein Binding:** {drug_info.get('protein_binding', 'N/A')}")
            st.markdown(f"- **Metabolism:** {drug_info.get('metabolism', 'N/A')}")
            
            if 'atc_codes' in drug_info and drug_info['atc_codes']:
                st.markdown("#### ATC Codes")
                for code in drug_info['atc_codes']:
                    st.code(code, language="text")

def main():
    """Main function to run the DrugBank viewer."""
    st.title("💊 DrugBank Viewer")
    
    # Example drug IDs
    example_drugs = {
        "Aspirin": "DB00945",
        "Ibuprofen": "DB01050",
        "Metformin": "DB00331",
        "Atorvastatin": "DB01076",
        "Omeprazole": "DB00338"
    }
    
    # Create two columns for the search interface
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Search box for drug ID
        drug_id = st.text_input(
            "Search DrugBank ID:", 
            value=st.query_params.get('drugbank', ''),
            placeholder="e.g., DB00173"
        )
    
    with col2:
        st.markdown("### ")
        if st.button("Search", use_container_width=True):
            if drug_id:
                st.query_params['drugbank'] = drug_id
            else:
                st.warning("Please enter a DrugBank ID")
    
    # Display example drugs
    st.markdown("### Example Drugs")
    cols = st.columns(len(example_drugs))
    for (name, db_id), col in zip(example_drugs.items(), cols):
        with col:
            if st.button(f"{name} ({db_id})", use_container_width=True, key=f"btn_{db_id}"):
                st.query_params['drugbank'] = db_id
                st.rerun()
    
    # Display drug information if ID is provided
    if 'drugbank' in st.query_params and st.query_params['drugbank']:
        drug_id = st.query_params['drugbank']
        display_drug_info(drug_id)
    else:
        st.info("🔍 Search for a drug using its DrugBank ID or click on an example above.")
        st.markdown("---")
        st.markdown("""
        ### How to use:
        1. Enter a DrugBank ID (e.g., DB00945 for Aspirin) in the search box
        2. Click the Search button or press Enter
        3. View detailed drug information and structure
        
        ### Example DrugBank IDs:
        - DB00945: Aspirin
        - DB01050: Ibuprofen
        - DB00331: Metformin
        - DB01076: Atorvastatin
        - DB00338: Omeprazole
        """)

if __name__ == "__main__":
    main()

# Try optional freesasa import (may be unavailable in some environments)
try:
    import freesasa
    FREESASA_AVAILABLE = True
except Exception:
    FREESASA_AVAILABLE = False

"""
PDB Processing and Visualization Utilities

This module provides functions for working with PDB files, including visualization,
ligand extraction, and structure analysis.
"""
import os
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import numpy as np
import pandas as pd
import py3Dmol
from stmol import showmol
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Atom import Atom
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_pdb(pdb_file: Union[str, Path], width: int = 800, height: int = 600) -> None:
    """
    Visualize a PDB file using py3Dmol.
    
    Args:
        pdb_file: Path to the PDB file
        width: Width of the visualization in pixels
        height: Height of the visualization in pixels
    """
    try:
        # Read PDB file
        with open(pdb_file, 'r') as f:
            pdb_data = f.read()
        
        # Create 3D visualization
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_data, 'pdb')
        
        # Style the protein and ligand
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.addStyle({'hetflag': True}, {'stick': {}})
        
        # Zoom to fit the structure
        view.zoomTo()
        view.zoom(1.2, 800)
        
        # Show in Streamlit
        showmol(view, height=height, width=width)
        
    except Exception as e:
        logger.error(f"Error visualizing PDB file: {str(e)}")
        raise

def extract_ligands(pdb_file: Union[str, Path]) -> List[Dict[str, any]]:
    """
    Extract ligand information from a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        List of dictionaries containing ligand information
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pdb', pdb_file)
        
        ligands = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if residue is a hetero residue (HETATM)
                    if residue.id[0].strip() != '':
                        # Get atoms for the ligand
                        atoms = []
                        for atom in residue:
                            atoms.append({
                                'name': atom.get_name(),
                                'element': atom.element,
                                'coords': atom.get_coord().tolist(),
                                'bfactor': atom.get_bfactor(),
                                'occupancy': atom.get_occupancy()
                            })
                        
                        # Add ligand info
                        ligands.append({
                            'residue_name': residue.get_resname(),
                            'chain_id': chain.id,
                            'residue_id': residue.id[1],
                            'insertion_code': residue.id[2],
                            'atoms': atoms,
                            'num_atoms': len(atoms)
                        })
        
        return ligands
    
    except Exception as e:
        logger.error(f"Error extracting ligands from PDB: {str(e)}")
        return []

def calculate_sasa(pdb_file: Union[str, Path]) -> Dict[str, float]:
    """
    Calculate Solvent Accessible Surface Area (SASA) for a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        Dictionary containing SASA values for the whole structure and by chain
    """
    try:
        if not FREESASA_AVAILABLE:
            logger.info("freesasa not available; returning default SASA values.")
            return {'total_sasa': 0.0, 'chain_sasa': {}}
        # Calculate SASA using FreeSASA
        structure = freesasa.Structure(str(pdb_file))
        result = freesasa.calc(structure)
        
        # Get total SASA
        total_sasa = result.totalArea()
        
        # Get SASA by chain
        chain_sasa = {}
        for chain_id in set(atom.chain for atom in structure.atoms()):
            chain_sasa[chain_id] = sum(
                atom.area for atom in structure.atoms() 
                if atom.chain == chain_id
            )
        
        return {
            'total_sasa': total_sasa,
            'chain_sasa': chain_sasa
        }
    
    except Exception as e:
        logger.error(f"Error calculating SASA: {str(e)}")
        return {'total_sasa': 0.0, 'chain_sasa': {}}

def clean_pdb(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    remove_hetatms: bool = False,
    keep_chains: List[str] = None,
    remove_water: bool = True
) -> bool:
    """
    Clean a PDB file by removing HETATMs, water molecules, and specified chains.
    
    Args:
        input_pdb: Path to the input PDB file
        output_pdb: Path to save the cleaned PDB file
        remove_hetatms: Whether to remove all HETATMs
        keep_chains: List of chain IDs to keep (None to keep all)
        remove_water: Whether to remove water molecules
        
    Returns:
        True if successful, False otherwise
    """
    class PDBSelect(Select):
        def accept_residue(self, residue):
            # Remove water if requested
            if remove_water and residue.get_resname() in ['HOH', 'WAT']:
                return 0
            
            # Remove HETATMs if requested
            if remove_hetatms and residue.id[0].strip() != '':
                return 0
                
            # Filter by chain if specified
            if keep_chains is not None and residue.get_parent().id not in keep_chains:
                return 0
                
            return 1
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pdb', input_pdb)
        
        # Clean the structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb), PDBSelect())
        
        return True
    
    except Exception as e:
        logger.error(f"Error cleaning PDB file: {str(e)}")
        return False

def extract_binding_pocket(
    pdb_file: Union[str, Path],
    ligand_residue: str,
    chain_id: str,
    distance: float = 5.0,
    output_file: Optional[Union[str, Path]] = None
) -> Optional[Path]:
    """
    Extract a binding pocket around a ligand.
    
    Args:
        pdb_file: Path to the PDB file
        ligand_residue: Residue ID of the ligand
        chain_id: Chain ID of the ligand
        distance: Distance in Angstroms around the ligand to include
        output_file: Path to save the extracted pocket (optional)
        
    Returns:
        Path to the extracted pocket PDB file or None if failed
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pdb', pdb_file)
        
        # Find the ligand
        ligand_atoms = []
        for model in structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                for residue in chain:
                    if str(residue.id[1]) == str(ligand_residue):
                        ligand_atoms.extend(atom for atom in residue)
        
        if not ligand_atoms:
            logger.error(f"Ligand {ligand_residue} not found in chain {chain_id}")
            return None
        
        # Get ligand coordinates
        ligand_coords = np.array([atom.get_coord() for atom in ligand_atoms])
        
        # Find residues within distance of the ligand
        pocket_residues = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Skip the ligand itself
                    if str(residue.id[1]) == str(ligand_residue) and chain.id == chain_id:
                        continue
                    
                    # Check distance to any ligand atom
                    for atom in residue:
                        if atom.element == 'H':  # Skip hydrogens
                            continue
                        
                        atom_coord = atom.get_coord()
                        dists = np.linalg.norm(ligand_coords - atom_coord, axis=1)
                        if np.any(dists <= distance):
                            pocket_residues.add((chain.id, residue.id[1]))
                            break
        
        # Create a new structure with just the pocket
        class PocketSelect(Select):
            def accept_residue(self, residue):
                return (residue.get_parent().id, residue.id[1]) in pocket_residues or \
                       (residue.get_parent().id == chain_id and str(residue.id[1]) == str(ligand_residue))
        
        # Save the pocket
        if output_file is None:
            output_file = Path(pdb_file).with_stem(f"{Path(pdb_file).stem}_pocket")
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_file), PocketSelect())
        
        return Path(output_file)
    
    except Exception as e:
        logger.error(f"Error extracting binding pocket: {str(e)}")
        return None

if __name__ == "__main__":
    # Example usage
    pdb_file = "example.pdb"
    
    # Visualize PDB
    visualize_pdb(pdb_file)
    
    # Extract ligands
    ligands = extract_ligands(pdb_file)
    print(f"Found {len(ligands)} ligands")
    
    # Calculate SASA
    sasa = calculate_sasa(pdb_file)
    print(f"Total SASA: {sasa['total_sasa']:.2f} Å²")
    
    # Clean PDB
    clean_pdb(pdb_file, "cleaned.pdb", remove_water=True)
    
    # Extract binding pocket
    if ligands:
        ligand = ligands[0]
        extract_binding_pocket(
            pdb_file,
            ligand_residue=ligand['residue_id'],
            chain_id=ligand['chain_id'],
            distance=5.0,
            output_file="pocket.pdb"
        )

"""Utility functions for the Drug Target Platform."""
import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import pandas as pd
import numpy as np
import plotly.express as px
from loguru import logger
import requests

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TargetProcessor:
    """Class for processing and analyzing target data."""
    
    def __init__(self, data_dir: Union[str, Path] = "data"):
        """Initialize with data directory."""
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
    
    def load_data(self, file_path: Union[str, Path, object]) -> pd.DataFrame:
        """
        Load data from CSV or Excel file.
        
        Args:
            file_path: Can be a string/Path to a file or a file-like object (e.g., from Streamlit's file_uploader)
            
        Returns:
            pd.DataFrame: Loaded data
        """
        try:
            # Handle file-like objects (e.g., from Streamlit's file_uploader)
            if hasattr(file_path, 'name') and hasattr(file_path, 'read'):
                file_name = file_path.name.lower()
                if file_name.endswith('.csv'):
                    return pd.read_csv(file_path)
                elif file_name.endswith(('.xlsx', '.xls')):
                    return pd.read_excel(file_path)
                else:
                    raise ValueError(f"Unsupported file format: {file_name.split('.')[-1]}")
            
            # Handle string/Path objects
            file_path = Path(file_path)
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")
                
            if file_path.suffix.lower() == '.csv':
                return pd.read_csv(file_path)
            elif file_path.suffix.lower() in ['.xlsx', '.xls']:
                return pd.read_excel(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_path.suffix}")
                
        except Exception as e:
            file_desc = getattr(file_path, 'name', str(file_path))
            logger.error(f"Error loading {file_desc}: {str(e)}")
            raise
    
    def preprocess_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Preprocess the input dataframe."""
        # Make a copy to avoid modifying the original
        df = df.copy()
        
        # Standardize column names (case-insensitive)
        col_map = {col: col.lower() for col in df.columns}
        df = df.rename(columns=col_map)
        
        # Required columns check
        required = {'gene', 'logfc', 'adj.p.val'}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        # Ensure numeric types
        df['logfc'] = pd.to_numeric(df['logfc'], errors='coerce')
        df['adj.p.val'] = pd.to_numeric(df['adj.p.val'], errors='coerce')
        
        # Remove rows with missing values in required columns
        df = df.dropna(subset=['gene', 'logfc', 'adj.p.val'])
        
        # Calculate score for prioritization
        df['score'] = abs(df['logfc']) * (-np.log10(df['adj.p.val']))
        
        return df
    
    def filter_targets(
        self, 
        df: pd.DataFrame, 
        min_logfc: float = 1.0, 
        max_pval: float = 0.05
    ) -> pd.DataFrame:
        """Filter targets based on logFC and p-value thresholds."""
        if df.empty:
            return df
            
        filtered = df[
            (abs(df['logfc']) >= min_logfc) & 
            (df['adj.p.val'] <= max_pval)
        ].copy()
        
        # Sort by score (descending)
        filtered = filtered.sort_values('score', ascending=False)
        
        return filtered
    
    def create_volcano_plot(
        self, 
        df: pd.DataFrame, 
        logfc_threshold: float = 1.0, 
        pval_threshold: float = 0.05
    ) -> px.scatter:
        """Create an interactive volcano plot."""
        if df.empty:
            return None
            
        # Create a copy to avoid modifying the original
        plot_df = df.copy()
        plot_df['-log10(pval)'] = -np.log10(plot_df['adj.p.val'])
        
        # Add significance status
        plot_df['significance'] = np.where(
            (abs(plot_df['logfc']) >= logfc_threshold) & 
            (plot_df['adj.p.val'] <= pval_threshold),
            'Significant',
            'Not Significant'
        )
        
        # Create the plot
        fig = px.scatter(
            plot_df, 
            x='logfc', 
            y='-log10(pval)',
            color='significance',
            color_discrete_map={
                'Significant': 'red',
                'Not Significant': 'grey'
            },
            hover_data=['gene'],
            labels={
                'logfc': 'log2 Fold Change',
                '-log10(pval)': '-log10(adj.P.Val)',
                'significance': 'Significance'
            },
            title='Volcano Plot of Differential Expression'
        )
        
        # Add threshold lines
        fig.add_hline(
            y=-np.log10(pval_threshold), 
            line_dash='dash', 
            line_color='black',
            annotation_text=f'p={pval_threshold}',
            annotation_position='bottom right'
        )
        
        for fc in [-logfc_threshold, logfc_threshold]:
            fig.add_vline(
                x=fc, 
                line_dash='dash', 
                line_color='black'
            )
        
        # Update layout
        fig.update_layout(
            legend_title_text='',
            showlegend=True,
            hovermode='closest'
        )
        
        return fig


def setup_directories(base_dir: Union[str, Path]) -> Dict[str, Path]:
    """Create and return a dictionary of directory paths."""
    base_dir = Path(base_dir)
    dirs = {
        'base': base_dir,
        'data': base_dir / 'data',
        'targets': base_dir / 'targets',
        'structures': base_dir / 'structures',
        'pockets': base_dir / 'pockets',
        'ligands': base_dir / 'ligands',
        'nanodelivery': base_dir / 'nanodelivery',
        'notebooks': base_dir / 'notebooks'
    }
    
    # Create directories if they don't exist
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    
    return dirs


def save_analysis_results(df: pd.DataFrame, filename: str, output_dir: Union[str, Path]) -> Path:
    """Save analysis results to a file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filepath = output_dir / filename
    
    if filename.endswith('.csv'):
        df.to_csv(filepath, index=False)
    elif filename.endswith(('.xlsx', '.xls')):
        df.to_excel(filepath, index=False)
    else:
        raise ValueError("Unsupported file format. Use .csv or .xlsx")
    
    logger.info(f"Results saved to {filepath}")
    return filepath


# ----------------------
# Protein sequence utils
# ----------------------
def fetch_uniprot_sequence(gene_symbol: str, organism_id: int = 9606) -> Optional[Tuple[str, str]]:
    """Fetch protein sequence from UniProt REST API for a given gene symbol.

    Args:
        gene_symbol: HGNC gene symbol (e.g., 'AMACR')
        organism_id: NCBI taxonomy ID (default: 9606 for human)

    Returns:
        (accession, sequence) if found, otherwise None
    """
    try:
        # Query UniProtKB for reviewed (Swiss-Prot) entries first, then fallback to all
        base = "https://rest.uniprot.org/uniprotkb/search"
        common_params = {
            "query": f"gene_exact:{gene_symbol} AND organism_id:{organism_id}",
            "fields": "accession,reviewed,protein_name,organism_id,length,sequence",
            "format": "json",
            "size": 1,
        }

        # Prefer reviewed entries
        params_reviewed = common_params.copy()
        params_reviewed["query"] += " AND reviewed:true"

        for params in (params_reviewed, common_params):
            resp = requests.get(base, params=params, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            results = data.get("results", [])
            if results:
                entry = results[0]
                accession = entry.get("primaryAccession")
                seq = entry.get("sequence", {}).get("value")
                if accession and seq:
                    return accession, seq
        return None
    except Exception as e:
        logger.error(f"UniProt fetch failed for {gene_symbol} ({organism_id}): {e}")
        return None


def save_fasta(name: str, sequence: str, out_path: Union[str, Path]) -> Path:
    """Save a sequence to FASTA file, ensuring parent directory exists.

    Args:
        name: header/name line to use
        sequence: amino acid sequence
        out_path: destination file path

    Returns:
        Path to written FASTA file
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">{name}\n")
        # Wrap sequence to 60 chars per line
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    logger.info(f"Saved FASTA to {out_path}")
    return out_path

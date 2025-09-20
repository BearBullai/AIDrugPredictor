# AI-Powered Drug Target Platform

A comprehensive platform for identifying and prioritizing drug targets, predicting protein structures, detecting binding pockets, and designing nanocarriers for targeted drug delivery.

## 🚀 Key Features

- **Data Processing**
  - Import and preprocess differential gene expression data
  - Interactive data visualization and filtering
  - Export results in multiple formats (CSV, Excel)

- **Target Analysis**
  - Volcano plots for differential expression analysis
  - Interactive filtering by log2 fold change and p-value
  - Detailed statistics and quality metrics

- **Structure Prediction**
  - High-quality 3D structure prediction using ColabFold/AlphaFold2
  - Visualize and analyze predicted structures
  - Confidence scoring with pLDDT metrics
  - Ligand binding site detection

- **Binding Pocket Analysis**
  - Automatic detection of potential binding pockets
  - Pocket druggability assessment
  - 3D visualization of binding sites

- **Ligand Docking**
  - Virtual screening of compound libraries
  - Molecular docking with AutoDock Vina
  - Binding affinity prediction
  - Pose visualization and analysis

- **Nanocarrier Design**
  - Lipid-based nanoparticle design
  - Polymer selection and optimization
  - Targeted delivery recommendations

- **Reporting & Export**
  - Generate comprehensive PDF reports
  - Export structures in PDB format
  - Download analysis results in multiple formats

## 🛠️ Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)
- Git (for cloning the repository)

### Setup Instructions

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd project
   ```

2. **Create and activate a virtual environment** (recommended):
   ```bash
   # Windows
   python -m venv venv
   .\venv\Scripts\activate
   
   # macOS/Linux
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
   
   Note: Some packages might require additional system dependencies. On Ubuntu/Debian, you might need:
   ```bash
   sudo apt-get install python3-dev python3-pip python3-venv
   ```

4. **Set up required directories**:
   ```bash
   mkdir -p data targets structures pockets ligands nanodelivery notebooks
   ```

## 🚀 Usage

1. **Start the application**:
   ```bash
   streamlit run app.py
   ```
   This will start the web server and open the application in your default browser at http://localhost:8501

2. **Upload your data**:
   - Use the "Data Upload" page to upload your differential gene expression data
   - Supported formats: CSV, Excel (XLSX), TSV

3. **Analyze targets**:
   - Filter and prioritize potential drug targets
   - View interactive visualizations of your data

4. **Predict structures**:
   - Select a target and predict its 3D structure
   - Analyze predicted structures and binding pockets

5. **Perform docking**:
   - Upload ligand structures or select from a library
   - Run molecular docking simulations
   - Analyze binding modes and interactions

6. **Design nanocarriers**:
   - Get recommendations for nanoparticle formulations
   - Optimize delivery parameters

## 📊 Example Workflow

1. Upload your differential gene expression data
2. Identify potential drug targets using the Target Analysis page
3. Predict 3D structures for top targets
4. Detect and analyze binding pockets
5. Perform virtual screening of compound libraries
6. Design nanocarriers for top hits
7. Generate comprehensive reports of your findings

## 📁 Project Structure

```
project/
├── data/                   # Input data files (DEG lists, etc.)
├── targets/                # Processed target information
├── structures/             # Predicted protein structures
│   ├── [target_name]/      # Per-target structure files
│   │   ├── [target].fasta  # Input sequence
│   │   ├── [target].pdb    # Predicted structure
│   │   └── [target]_plddt.json  # Confidence scores
├── pockets/                # Detected binding pockets
├── ligands/                # Ligand structures and docking results
├── nanodelivery/           # Nanocarrier design files
├── notebooks/              # Jupyter notebooks for analysis
├── pages/                  # Streamlit page modules
│   ├── 1_Data_Upload.py    # Data import and preprocessing
│   ├── 2_Target_Analysis.py # Target prioritization
│   ├── 3_Structure_Prediction.py # 3D structure prediction
│   ├── 4_Pocket_Detection.py # Binding pocket analysis
│   └── ...
├── utils/                  # Utility functions
├── config.py              # Configuration settings
├── pdb_utils.py           # PDB processing utilities
└── requirements.txt       # Python dependencies
├── ligands/                # Ligand data and docking results
├── nanodelivery/           # Nanocarrier design specifications
├── notebooks/              # Jupyter notebooks for analysis
├── app.py                  # Main Streamlit application
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

## Data Requirements

The input data file should contain at least the following columns:
- `gene`: Gene symbol
- `logFC`: Log2 fold change
- `adj.P.Val`: Adjusted p-value

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Streamlit for the web framework
- ColabFold/AlphaFold2 for structure prediction
- RDKit for cheminformatics
- Plotly for interactive visualizations

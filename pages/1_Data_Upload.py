"""Data Upload Page"""
import streamlit as st
import pandas as pd
from pathlib import Path
import logging
from utils import TargetProcessor
from config import config, ensure_dir

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set page config
st.set_page_config(
    page_title="Data Upload - AI Drug Target Platform",
    page_icon="📤"
)

# Initialize session state
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'df' not in st.session_state:
    st.session_state.df = None
if 'filtered_df' not in st.session_state:
    st.session_state.filtered_df = None
if 'processor' not in st.session_state:
    st.session_state.processor = TargetProcessor(Path("data"))

# Main function
def main():
    st.title("📤 Data Upload")
    st.markdown("---")
    
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
            
            # Show data preview with more rows and better formatting
            st.subheader("📊 Data Overview")
            
            # Basic stats in a more comprehensive way
            stats_col1, stats_col2, stats_col3, stats_col4 = st.columns(4)
            with stats_col1:
                st.metric("Total Genes", len(df))
            with stats_col2:
                st.metric("Upregulated", f"{sum(df['logfc'] > 0):,}")
            with stats_col3:
                st.metric("Downregulated", f"{sum(df['logfc'] < 0):,}")
            with stats_col4:
                st.metric("Significant (p < 0.05)", f"{sum(df['adj.p.val'] < 0.05):,}")
            
            # Add distribution plots
            st.subheader("📈 Data Distribution")
            plot_col1, plot_col2 = st.columns(2)
            
            with plot_col1:
                st.bar_chart(data=df['logfc'].value_counts().sort_index(), 
                           use_container_width=True)
                st.caption("LogFC Distribution")
                
            with plot_col2:
                st.line_chart(data=df['adj.p.val'].value_counts().sort_index(), 
                            use_container_width=True)
                st.caption("Adjusted p-value Distribution")
            
            # Enhanced data preview with search and sort
            st.subheader("🔍 Data Preview")
            
            # Add search functionality
            search_term = st.text_input("Search in data", "", 
                                      placeholder="Type to filter rows...",
                                      help="Search across all columns")
            
            # Filter data based on search
            if search_term:
                try:
                    # Try to convert search term to number for numeric comparison
                    try:
                        search_num = float(search_term)
                        filtered_df = df[df.apply(lambda row: row.astype(str).str.contains(search_term, case=False).any() | 
                                                (row == search_num).any(), axis=1)]
                    except ValueError:
                        filtered_df = df[df.apply(lambda row: row.astype(str).str.contains(search_term, case=False).any(), axis=1)]
                except Exception as e:
                    st.warning(f"Error in search: {str(e)}")
                    filtered_df = df
            else:
                filtered_df = df
            
            # Display the first 100 rows with sorting
            st.dataframe(
                filtered_df.head(100),
                use_container_width=True,
                height=400,
                hide_index=True,
                column_config={
                    "logfc": st.column_config.NumberColumn(
                        "LogFC",
                        format="%.3f",
                        help="Log2 Fold Change"
                    ),
                    "adj.p.val": st.column_config.NumberColumn(
                        "Adjusted p-value",
                        format="%.3e",
                        help="Adjusted p-value"
                    )
                }
            )
            
            # Show data summary
            with st.expander("📋 Data Summary", expanded=False):
                st.write(filtered_df.describe())
                
            # Download button for filtered data
            if len(filtered_df) < len(df):
                csv = filtered_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="💾 Download Filtered Data",
                    data=csv,
                    file_name="filtered_deg_data.csv",
                    mime="text/csv"
                )
            
            # Save the uploaded file
            save_path = Path("data") / "deg_list.csv"
            save_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(save_path, index=False)
            
            # Show file info
            st.info(f"💾 Data saved to: `{save_path}`")
            
        except Exception as e:
            st.error(f"❌ Error processing file: {str(e)}")
            logger.exception("Error in data upload")
    else:
        st.info("ℹ️ Please upload a file to get started")
        
        # Check if data already exists
        default_data = Path("data") / "deg_list.csv"
        if default_data.exists():
            if st.button("📂 Load existing data"):
                try:
                    df = pd.read_csv(default_data)
                    df = st.session_state.processor.preprocess_data(df)
                    
                    st.session_state.df = df
                    st.session_state.data_loaded = True
                    st.session_state.filtered_df = None
                    
                    st.success("✅ Existing data loaded successfully!")
                    st.dataframe(df.head())
                except Exception as e:
                    st.error(f"❌ Error loading existing data: {str(e)}")
                    logger.exception("Error loading existing data")

# Run the app
if __name__ == "__main__":
    main()

"""Target Analysis Page"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import logging
import io
from utils import TargetProcessor
from config import config
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set page config
st.set_page_config(
    page_title="Target Analysis - AI Drug Target Platform",
    page_icon="🔍"
)

# Check if data is loaded
if 'data_loaded' not in st.session_state or not st.session_state.data_loaded:
    st.warning("⚠️ Please upload data in the 'Data Upload' section first.")
    if st.button("Go to Data Upload"):
        st.switch_page("pages/1_Data_Upload.py")
    st.stop()

# Main function
def main():
    st.title("🔍 Target Analysis")
    st.markdown("---")
    
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
        
        # Add DrugBank compound accession codes if not present
        if 'drugbank_id' not in filtered_df.columns:
            # This is a placeholder - in a real app, you would fetch this from a database
            filtered_df['drugbank_id'] = 'DB' + (filtered_df.index + 10000).astype(str)
        
        # Show top targets with more columns and better formatting
        display_cols = ['gene', 'logfc', 'adj.p.val', 'score', 'drugbank_id']
        display_cols = [col for col in display_cols if col in filtered_df.columns]
        
        # Display the dataframe with better formatting
        st.data_editor(
            filtered_df[display_cols].head(top_n),
            column_config={
                'gene': st.column_config.TextColumn('Gene', width='small'),
                'logfc': st.column_config.NumberColumn(
                    'Log2 FC',
                    format='%.2f',
                    help='Log2 Fold Change',
                    width='small'
                ),
                'adj.p.val': st.column_config.NumberColumn(
                    'Adj. p-value',
                    format='%.3e',
                    help='Adjusted p-value',
                    width='small'
                ),
                'score': st.column_config.ProgressColumn(
                    'Score',
                    help='Target score',
                    format='%.2f',
                    min_value=0,
                    max_value=float(filtered_df['score'].max()) * 1.1,
                    width='medium'
                ),
                'drugbank_id': st.column_config.LinkColumn(
                    'DrugBank',
                    help='DrugBank compound accession',
                    display_text='View',
                    width='small'
                )
            },
            hide_index=True,
            use_container_width=True,
            height=min(400, 35 * min(len(filtered_df), 20)),
            column_order=display_cols
        )
        
        # Add pagination for large datasets
        if len(filtered_df) > top_n:
            st.caption(f"Showing top {top_n} of {len(filtered_df)} targets. Use filters to narrow down results.")
        
        # Download buttons
        col1, col2 = st.columns(2)
        with col1:
            csv = filtered_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Full Results (CSV)",
                data=csv,
                file_name="filtered_targets.csv",
                mime="text/csv",
                use_container_width=True
            )
        with col2:
            excel = filtered_df.to_excel("filtered_targets.xlsx", index=False, engine='openpyxl')
            with open("filtered_targets.xlsx", "rb") as f:
                st.download_button(
                    label="📥 Download Full Results (Excel)",
                    data=f,
                    file_name="filtered_targets.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
        
        # Show detailed statistics
        st.subheader("📈 Detailed Statistics")
        
        # Create tabs for different statistics
        tab1, tab2, tab3 = st.tabs(["Summary", "Distribution", "Top Targets"])
        
        with tab1:
            # Basic statistics
            stats_cols = st.columns(3)
            with stats_cols[0]:
                st.metric("Total Targets", len(filtered_df))
                st.metric("Up-regulated", len(filtered_df[filtered_df['logfc'] > 0]))
            with stats_cols[1]:
                st.metric("Mean Score", f"{filtered_df['score'].mean():.2f}")
                st.metric("Down-regulated", len(filtered_df[filtered_df['logfc'] < 0]))
            with stats_cols[2]:
                st.metric("Top Score", f"{filtered_df['score'].max():.2f}")
                st.metric("Avg. |Log2FC|", f"{filtered_df['logfc'].abs().mean():.2f}")
        
        with tab2:
            # Distribution plots
            col1, col2 = st.columns(2)
            with col1:
                st.plotly_chart(
                    px.histogram(
                        filtered_df, 
                        x='logfc',
                        title='Log2 Fold Change Distribution',
                        labels={'logfc': 'Log2 Fold Change'},
                        color_discrete_sequence=['#636EFA']
                    ),
                    use_container_width=True
                )
            with col2:
                st.plotly_chart(
                    px.histogram(
                        filtered_df, 
                        x='score',
                        title='Score Distribution',
                        labels={'score': 'Target Score'},
                        color_discrete_sequence=['#EF553B']
                    ),
                    use_container_width=True
                )
        
        with tab3:
            # Top targets table with more details
            top_targets = filtered_df.nlargest(10, 'score')
            st.dataframe(
                top_targets[['gene', 'logfc', 'adj.p.val', 'score', 'drugbank_id']],
                column_config={
                    'gene': 'Gene',
                    'logfc': st.column_config.NumberColumn('Log2 FC', format='%.2f'),
                    'adj.p.val': st.column_config.NumberColumn('Adj. p-value', format='%.2e'),
                    'score': st.column_config.ProgressColumn('Score', format='%.2f', min_value=0, max_value=1),
                    'drugbank_id': 'DrugBank ID'
                },
                use_container_width=True,
                hide_index=True
            )
    else:
        st.info("ℹ️ Apply filters to see the filtered results")

# Run the app
if __name__ == "__main__":
    main()

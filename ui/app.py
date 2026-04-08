import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import glob

# Page configuration
st.set_page_config(page_title="NumPhys: Transient Results Visualizer", layout="wide")

st.title("📂 NumPhys Simulation Visualizer")
st.markdown("""
This dashboard allows you to explore **transient simulation results** stored in the `backend/output` directory.
Select a result set below to visualize how the physics evolved across both **space and time**.
""")

# Output directory path
OUTPUT_DIR = os.path.abspath(os.path.join(os.getcwd(), "backend", "output"))

# Helper to find result files
def get_result_files():
    if not os.path.exists(OUTPUT_DIR):
        return []
    return glob.glob(os.path.join(OUTPUT_DIR, "*.csv"))

files = get_result_files()

if not files:
    st.warning(f"No simulation results found in `{OUTPUT_DIR}`. Please run the backend simulations first!")
    st.stop()

# Sidebar: File Selection
selected_file = st.sidebar.selectbox(
    "Select Simulation Result",
    files,
    format_func=lambda x: os.path.basename(x)
)

# Load data
@st.cache_data
def load_data(file_path):
    df = pd.read_csv(file_path)
    return df

try:
    df = load_data(selected_file)
    # Detect value column (temperature, pressure, etc.)
    val_col = [col for col in df.columns if col not in ['time', 'index', 'x']][0]
    
    times = sorted(df['time'].unique())
    
    st.sidebar.divider()
    st.sidebar.header("Animation Controls")
    
    # Time scrubbing
    selected_time = st.sidebar.select_slider(
        "Current Time",
        options=times,
        format_func=lambda x: f"t = {x:.2f}"
    )

    # Plot spatial profile at selected time
    sub_df = df[df['time'] == selected_time]
    
    st.subheader(f"📊 Spatial Distribution ({val_col.capitalize()})")
    
    fig = px.line(
        sub_df, 
        x='x', 
        y=val_col, 
        title=f"{val_col.capitalize()} Profile at t = {selected_time:.2f}",
        labels={'x': 'Position (x)', val_col: val_col.capitalize()},
        markers=True,
        template="plotly_dark"
    )
    
    fig.update_layout(
        yaxis_range=[df[val_col].min() * 0.9, df[val_col].max() * 1.1],
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    st.plotly_chart(fig, use_container_width=True)

    # Secondary Plot: Time evolution at a specific point
    st.divider()
    st.subheader("⏱️ Temporal Evolution")
    
    # Select point to watch
    indices = sorted(df['index'].unique())
    selected_idx = st.selectbox("Select grid point to monitor", options=indices, format_func=lambda x: f"Node {x} (x={df[df['index']==x]['x'].iloc[0]:.2f})")
    
    point_data = df[df['index'] == selected_idx]
    
    fig_time = px.line(
        point_data,
        x='time',
        y=val_col,
        title=f"{val_col.capitalize()} over Time at Node {selected_idx}",
        template="plotly_dark",
        labels={'time': 'Time (t)', val_col: val_col.capitalize()}
    )
    st.plotly_chart(fig_time, use_container_width=True)

    # Animated Heatmap (Space vs Time)
    st.divider()
    st.subheader("🌊 Space-Time Heatmap")
    
    # Pivot for heatmap
    pivot_df = df.pivot(index='time', columns='x', values=val_col)
    
    fig_heat = px.imshow(
        pivot_df,
        labels=dict(x="Position (x)", y="Time (t)", color=val_col.capitalize()),
        x=pivot_df.columns,
        y=pivot_df.index,
        color_continuous_scale='Viridis',
        aspect="auto",
        template="plotly_dark"
    )
    st.plotly_chart(fig_heat, use_container_width=True)

except Exception as e:
    st.error(f"Error loading or processing file: {e}")

st.divider()
st.caption("Axscnt NumPhys Project - Backend Results Visualizer")

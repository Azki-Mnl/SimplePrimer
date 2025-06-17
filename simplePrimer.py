import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from io import StringIO
import pandas as pd
import plotly.graph_objects as go
import re
import time
import os
from itertools import product
from collections import defaultdict

# --- Constants ---
DEFAULT_PARAMS = {
    'primer_len': 20,
    'gc_min': 40,
    'gc_max': 60,
    'tm_min': 50,
    'tm_max': 60,
    'min_amplicon': 100,
    'max_amplicon': 1000,
    'max_dimer_score': 4,
    'min_primer_distance': 50
}

EXAMPLE_SEQUENCE = """>Example DNA Sequence
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
"""

# --- Helper Functions ---
def is_valid_dna(seq):
    """Check if sequence contains only valid DNA characters"""
    if not seq:
        return False
    clean_seq = ''.join([line.strip().upper() for line in seq.splitlines() if not line.startswith(">")])
    return bool(re.fullmatch(r"^[ATCG]+$", clean_seq))

def reverse_complement(seq):
    """Generate reverse complement of DNA sequence"""
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def calculate_dimer_score(seq1, seq2):
    """Calculate potential dimerization score between two sequences"""
    max_score = 0
    for i in range(1, min(len(seq1), len(seq2)) + 1):
        # Check end of seq1 vs start of seq2
        score = sum(a == b for a, b in zip(seq1[-i:], seq2[:i]))
        max_score = max(max_score, score)
        
        # Check end of seq2 vs start of seq1
        score = sum(a == b for a, b in zip(seq2[-i:], seq1[:i]))
        max_score = max(max_score, score)
    return max_score

def find_primers(seq, primer_len, gc_min, gc_max, tm_min, tm_max):
    """Find all potential primers in sequence meeting criteria"""
    results = []
    for i in range(len(seq) - primer_len + 1):
        primer = seq[i:i+primer_len]
        gc = gc_fraction(primer) * 100
        tm = mt.Tm_Wallace(primer)
        if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
            results.append({
                'pos': i,
                'seq': primer,
                'gc': gc,
                'tm': tm,
                'end': i + primer_len
            })
    return results

def find_optimal_pairs(fwd_primers, rev_primers, seq_len, params):
    """Find optimal primer pairs with sequential region selection and best pairs first"""
    # First find all possible valid pairs with their scores
    all_pairs = []
    for fwd in fwd_primers:
        for rev in rev_primers:
            amp_size = (seq_len - rev['pos'] - params['primer_len']) - fwd['pos']
            
            # Check basic constraints
            if not (params['min_amplicon'] <= amp_size <= params['max_amplicon']):
                continue
                
            # Check dimer formation
            dimer_score = calculate_dimer_score(fwd['seq'], rev['seq'])
            if dimer_score > params['max_dimer_score']:
                continue
                
            # Calculate pair properties
            tm_diff = abs(fwd['tm'] - rev['tm'])
            gc_diff = abs(fwd['gc'] - rev['gc'])
            
            # Calculate composite score (lower is better)
            score = (dimer_score * 0.4) + (tm_diff * 0.3) + (gc_diff * 0.3)
            
            all_pairs.append({
                'fwd': fwd,
                'rev': rev,
                'amp_size': amp_size,
                'dimer_score': dimer_score,
                'tm_diff': tm_diff,
                'gc_diff': gc_diff,
                'score': score,
                'fwd_pos': fwd['pos'],
                'rev_pos': seq_len - rev['end']
            })
    
    # Sort all pairs by score (best first)
    all_pairs.sort(key=lambda x: x['score'])
    
    # Select best non-overlapping pairs
    selected_pairs = []
    used_positions = []
    min_distance = params['min_primer_distance']
    
    for pair in all_pairs:
        # Check if this pair overlaps with any used regions
        overlap = False
        fwd_start = pair['fwd_pos']
        fwd_end = fwd_start + params['primer_len']
        rev_start = pair['rev_pos']
        rev_end = rev_start + params['primer_len']
        
        for used_start, used_end in used_positions:
            if (fwd_start < used_end + min_distance and fwd_end > used_start - min_distance) or \
               (rev_start < used_end + min_distance and rev_end > used_start - min_distance):
                overlap = True
                break
                
        if not overlap:
            selected_pairs.append(pair)
            used_positions.append((fwd_start, fwd_end))
            used_positions.append((rev_start, rev_end))
            
            if len(selected_pairs) >= 10:  # Limit to top 10 pairs
                break
    
    # Resort the selected pairs by position while maintaining quality
    selected_pairs.sort(key=lambda x: x['fwd_pos'])
    
    return selected_pairs

def style_dataframe(df):
    """Style the results dataframe with highlights and gradients"""
    # Find best scores
    best_quality = df['Quality Score'].min()
    best_dimer = df['Dimer Score'].min()
    
    def highlight_row(row):
        styles = [''] * len(row)
        if row['Quality Score'] == best_quality:
            styles[df.columns.get_loc('Quality Score')] = 'background-color: yellow'
        if row['Dimer Score'] == best_dimer:
            styles[df.columns.get_loc('Dimer Score')] = 'background-color: lightgreen'
        return styles
    
    styled_df = df.style.apply(highlight_row, axis=1)
    
    # Add gradients (more compatible version)
    try:
        styled_df = styled_df.background_gradient(
            subset=['Quality Score'], 
            cmap='YlOrRd_r',
            vmin=df['Quality Score'].min(),
            vmax=df['Quality Score'].max()
        )
        styled_df = styled_df.background_gradient(
            subset=['Dimer Score'], 
            cmap='Greens_r',
            vmin=df['Dimer Score'].min(),
            vmax=df['Dimer Score'].max()
        )
    except:
        # Fallback if gradient fails
        pass
    
    return styled_df

def plot_primers(seq_len, pairs, primer_len):
    """Visualize primer binding positions"""
    fig = go.Figure()
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for idx, pair in enumerate(pairs[:5]):  # Show top 5 pairs
        color = colors[idx % len(colors)]
        pair_num = f"Pair {idx+1}"
        
        # Forward primer
        fig.add_trace(go.Scatter(
            x=[pair['fwd']['pos'], pair['fwd']['end']],
            y=[idx, idx],
            mode='lines+text',
            line=dict(color=color, width=8),
            text=[None, f"F{idx+1}"],
            textposition="top center",
            hoverinfo="text",
            hovertext=(
                f"{pair_num} - Forward<br>"
                f"Pos: {pair['fwd']['pos']}-{pair['fwd']['end']}<br>"
                f"Seq: {pair['fwd']['seq']}<br>"
                f"GC: {pair['fwd']['gc']:.1f}%<br>"
                f"Tm: {pair['fwd']['tm']:.1f}°C"
            )
        ))
        
        # Reverse primer (position from end)
        rev_pos = seq_len - pair['rev']['end']
        rev_end = seq_len - pair['rev']['pos']
        fig.add_trace(go.Scatter(
            x=[rev_pos, rev_end],
            y=[idx, idx],
            mode='lines+text',
            line=dict(color=color, width=8, dash='dot'),
            text=[None, f"R{idx+1}"],
            textposition="top center",
            hoverinfo="text",
            hovertext=(
                f"{pair_num} - Reverse<br>"
                f"Pos: {rev_pos}-{rev_end}<br>"
                f"Seq: {pair['rev']['seq']}<br>"
                f"GC: {pair['rev']['gc']:.1f}%<br>"
                f"Tm: {pair['rev']['tm']:.1f}°C"
            )
        ))
        
        # Amplicon line
        fig.add_trace(go.Scatter(
            x=[pair['fwd']['end'], rev_pos],
            y=[idx, idx],
            mode='lines',
            line=dict(color=color, width=1, dash='dash'),
            hoverinfo="none",
            showlegend=False
        ))
    
    fig.update_layout(
        height=300 + min(len(pairs), 5) * 30,
        title="Primer Binding Positions and Amplicons",
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(showticklabels=False),
        showlegend=False,
        hovermode="closest",
        margin=dict(l=20, r=20, t=60, b=20),
        template="plotly_white"
    )
    return fig

def reset_session():
    """Reset all session state variables"""
    st.session_state.clear()
    st.session_state.update({
        'seq': None,
        'primers_ready': False,
        'pairs': None,
        'params': DEFAULT_PARAMS.copy(),
        'widget_keys': {
            'text_area': f"text_area_{time.time()}",
            'file_uploader': f"file_uploader_{time.time()}"
        },
        'seq_auto_submitted': False
    })

# --- UI Configuration ---
def setup_page():
    """Configure page settings and styles"""
    st.set_page_config(
        page_title="SimplePrimer",
        layout="wide",
        page_icon="🧬" 
    )
    
    # Try to load background image
    bg_base64 = ""
    if os.path.exists("bg.txt"):
        try:
            with open("bg.txt", "r") as f:
                bg_base64 = f.read()
        except Exception:
            pass
    
    # Custom styling
    st.markdown(f"""
        <style>
        .stApp {{
            background-image: url("data:image/png;base64,{bg_base64}");
            background-size: cover;
            background-position: center;
            background-repeat: no-repeat;
            background-attachment: fixed;
        }}
        .css-18e3th9 {{
            background-color: rgba(255, 255, 255, 0.9);
            border-radius: 10px;
            padding: 2rem;
        }}
        .st-bb, .st-at, .st-ae, .st-af, .st-ag, .st-ah, .st-ai, .st-aj {{
            background-color: rgba(255, 255, 255, 0.8) !important;
        }}
        .stTextArea textarea, .stTextInput input {{
            background-color: rgba(240, 255, 255, 0.9) !important;
        }}
        .stDownloadButton, .stButton>button {{
            background-color: #140b42 !important;
            color: Blue !important;
            border-radius: 5px;
            border: 1px solid #444;
        }}
        .stAlert {{
            border-radius: 5px;
        }}
        /* Always show button text with gold color */
        .stDownloadButton > button span,
        .stButton > button span {{
            visibility: visible !important;
            color: #FFD700 !important;  /* Gold color */
            font-weight: 500;
        }}
        /* Button hover effects */
        .stDownloadButton > button:hover,
        .stButton > button:hover {{
            background-color: #2a1a6e !important;
            border-color: #FFD700;
        }}
        /* Active/clicked button effect */
        .stDownloadButton > button:active,
        .stButton > button:active {{
            background-color: #1a1042 !important;
        }}
        </style>
    """, unsafe_allow_html=True)
    
    # App title
    st.markdown("""
        <h1 style="text-align:center; color: #080740; margin-bottom: 0.5rem;">SimplePrimer</h1>
        <p style="text-align:center; font-size: 1.2rem; color: #13131f; margin-bottom: 2rem;">
        Automated PCR Primer Design</p>
    """, unsafe_allow_html=True)

# --- Main App ---
def main():
    setup_page()
    
    # Initialize session state
    if 'params' not in st.session_state:
        st.session_state.params = DEFAULT_PARAMS.copy()
    
    if 'widget_keys' not in st.session_state:
        st.session_state.widget_keys = {
            'text_area': f"text_area_{time.time()}",
            'file_uploader': f"file_uploader_{time.time()}"
        }
    
    if 'seq_auto_submitted' not in st.session_state:
        st.session_state.seq_auto_submitted = False
    
    # Layout columns
    left, right = st.columns([1, 1.2], gap="large")
    
    with left:
        st.subheader("📋 Sequence Input")
        
        # Example sequence button
        if st.button("💡 Load Example Sequence", help="Load an example DNA sequence to test the tool"):
            st.session_state.widget_keys = {
                'text_area': f"text_area_{time.time()}",
                'file_uploader': f"file_uploader_{time.time()}"
            }
            st.session_state.seq = ''.join(
                [line.strip().upper() for line in EXAMPLE_SEQUENCE.splitlines() 
                 if not line.startswith(">")]
            )
            st.session_state.primers_ready = False
            st.session_state.seq_auto_submitted = True
            st.rerun()
        
        # Sequence input with managed key
        fasta_input = st.text_area(
            label="Paste your FASTA sequence here:",
            height=200,
            key=st.session_state.widget_keys['text_area'],
            help="DNA sequence in FASTA format (only A, T, C, G characters allowed)"
        )
        
        # Add submit button
        submit_seq = st.button("Submit Sequence", type="primary")
        
        # Input validation (modified to check both the text area and submit button)
        if fasta_input and (submit_seq or st.session_state.get('seq_auto_submitted')):
            if is_valid_dna(fasta_input):
                clean_seq = ''.join([line.strip().upper() for line in fasta_input.splitlines() 
                                   if not line.startswith(">")])
                if st.session_state.get('seq') != clean_seq:
                    st.session_state.seq = clean_seq
                    st.session_state.primers_ready = False
                    st.success("✅ Valid DNA sequence loaded")
            else:
                st.error("❌ Invalid DNA sequence. Only A, T, C, G characters allowed.")
        
        # File uploader with managed key
        uploaded_file = st.file_uploader(
            "Or upload a FASTA file",
            type=["fasta", "fa", "txt"],
            key=st.session_state.widget_keys['file_uploader'],
            help="Upload a FASTA format file containing your DNA sequence"
        )
        
        if uploaded_file:
            try:
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                record = SeqIO.read(stringio, "fasta")
                seq = str(record.seq).upper()
                if is_valid_dna(seq):
                    st.session_state.seq = seq
                    st.session_state.primers_ready = False
                    st.session_state.seq_auto_submitted = True
                    st.success("✅ Valid DNA sequence loaded from file")
                else:
                    st.error("❌ Uploaded file contains invalid DNA sequence")
            except Exception as e:
                st.error(f"❌ Error reading file: {str(e)}")
        
        # Clear input button
        if st.button("❌ Clear Input", help="Reset all inputs and results"):
            reset_session()
            st.rerun()
    
    with right:
        st.subheader("🧪 Primer Design")
        
        if st.session_state.get('seq'):
            # Display sequence info
            seq_len = len(st.session_state.seq)
            st.info(f"Loaded sequence: {seq_len} bp | GC: {gc_fraction(st.session_state.seq)*100:.1f}%")
            
            # Mode selection
            mode = st.radio(
                "Design Mode:",
                ["Quick Design (Recommended)", "Advanced Parameters"],
                horizontal=True
            )
            
            # Parameter controls
            if mode == "Advanced Parameters":
                with st.expander("⚙️ Primer Design Parameters", expanded=True):
                    cols = st.columns(2)
                    with cols[0]:
                        st.session_state.params['primer_len'] = st.slider(
                            "Primer Length (bp)",
                            16, 30, DEFAULT_PARAMS['primer_len'],
                            help="Length of primers to design (typically 18-22 bp)"
                        )
                        st.session_state.params['gc_min'] = st.slider(
                            "Minimum GC%",
                            30, 60, DEFAULT_PARAMS['gc_min'],
                            help="Minimum GC content percentage"
                        )
                        st.session_state.params['tm_min'] = st.slider(
                            "Minimum Tm (°C)",
                            45, 60, DEFAULT_PARAMS['tm_min'],
                            help="Minimum melting temperature"
                        )
                        st.session_state.params['min_amplicon'] = st.slider(
                            "Minimum Amplicon Size (bp)",
                            50, 500, DEFAULT_PARAMS['min_amplicon'],
                            help="Minimum PCR product size"
                        )
                    
                    with cols[1]:
                        st.session_state.params['gc_max'] = st.slider(
                            "Maximum GC%",
                            st.session_state.params['gc_min']+1, 80, DEFAULT_PARAMS['gc_max'],
                            help="Maximum GC content percentage"
                        )
                        st.session_state.params['tm_max'] = st.slider(
                            "Maximum Tm (°C)",
                            st.session_state.params['tm_min']+1, 80, DEFAULT_PARAMS['tm_max'],
                            help="Maximum melting temperature"
                        )
                        st.session_state.params['max_amplicon'] = st.slider(
                            "Maximum Amplicon Size (bp)",
                            st.session_state.params['min_amplicon']+50, 2000, DEFAULT_PARAMS['max_amplicon'],
                            help="Maximum PCR product size"
                        )
                        st.session_state.params['max_dimer_score'] = st.slider(
                            "Maximum Dimer Score",
                            0, 8, DEFAULT_PARAMS['max_dimer_score'],
                            help="Lower values reduce primer-dimer potential"
                        )
                        st.session_state.params['min_primer_distance'] = st.slider(
                            "Minimum Primer Distance (bp)",
                            20, 200, DEFAULT_PARAMS['min_primer_distance'],
                            help="Minimum distance between primer pairs"
                        )
            
            # Design primers button
            if st.button("🧬 Design Primers", type="primary", help="Find optimal primer pairs"):
                if st.session_state.seq:
                    with st.spinner("Searching for optimal primer pairs..."):
                        try:
                            # Find forward and reverse primers
                            fwd_primers = find_primers(
                                st.session_state.seq,
                                st.session_state.params['primer_len'],
                                st.session_state.params['gc_min'],
                                st.session_state.params['gc_max'],
                                st.session_state.params['tm_min'],
                                st.session_state.params['tm_max']
                            )
                            
                            rev_primers = find_primers(
                                reverse_complement(st.session_state.seq),
                                st.session_state.params['primer_len'],
                                st.session_state.params['gc_min'],
                                st.session_state.params['gc_max'],
                                st.session_state.params['tm_min'],
                                st.session_state.params['tm_max']
                            )
                            
                            # Find optimal pairs
                            pairs = find_optimal_pairs(
                                fwd_primers,
                                rev_primers,
                                len(st.session_state.seq),
                                st.session_state.params
                            )
                            
                            if pairs:
                                st.session_state.pairs = pairs[:10]  # Keep top 10 pairs
                                st.session_state.primers_ready = True
                                st.success(f"Found {len(pairs)} potential primer pairs")
                            else:
                                st.session_state.primers_ready = False
                                st.error("No suitable primer pairs found. Try adjusting parameters.")
                        
                        except Exception as e:
                            st.error(f"Error during primer design: {str(e)}")
                            st.session_state.primers_ready = False
            
            # Display results if available
            if st.session_state.get('primers_ready', False) and st.session_state.get('pairs'):
                st.subheader("🔍 Results (Best Primers First)")
                
                # Create results table
                results = []
                for i, pair in enumerate(st.session_state.pairs, 1):
                    results.append({
                        'Pair': f"Pair {i}",
                        'Fwd Sequence': pair['fwd']['seq'],
                        'Fwd Position': f"{pair['fwd']['pos']}-{pair['fwd']['end']}",
                        'Fwd GC%': f"{pair['fwd']['gc']:.1f}",
                        'Fwd Tm (°C)': f"{pair['fwd']['tm']:.1f}",
                        'Rev Sequence': pair['rev']['seq'],
                        'Rev Position': f"{len(st.session_state.seq)-pair['rev']['end']}-{len(st.session_state.seq)-pair['rev']['pos']}",
                        'Rev GC%': f"{pair['rev']['gc']:.1f}",
                        'Rev Tm (°C)': f"{pair['rev']['tm']:.1f}",
                        'Amplicon Size': pair['amp_size'],
                        'Dimer Score': pair['dimer_score'],
                        'Quality Score': f"{pair['score']:.2f}"
                    })
                
                df = pd.DataFrame(results)
                
                # Display the styled dataframe
                try:
                    st.dataframe(
                        style_dataframe(df),
                        use_container_width=True
                    )
                except Exception as e:
                    st.warning("Couldn't apply advanced styling to the table")
                    st.dataframe(df)  # Fallback to basic table
                
                # Create FASTA output
                fasta_out = []
                for i, pair in enumerate(st.session_state.pairs, 1):
                    fasta_out.append(f">Pair{i}_Fwd_Pos{pair['fwd']['pos']}-{pair['fwd']['end']}_Score{pair['score']:.2f}")
                    fasta_out.append(pair['fwd']['seq'])
                    fasta_out.append(f">Pair{i}_Rev_Pos{len(st.session_state.seq)-pair['rev']['end']}-{len(st.session_state.seq)-pair['rev']['pos']}_Score{pair['score']:.2f}")
                    fasta_out.append(pair['rev']['seq'])
                
                # Download buttons
                col1, col2 = st.columns(2)
                with col1:
                    st.download_button(
                        "📥 Download Primer Pairs (CSV)",
                        df.to_csv(index=False),
                        "primers.csv",
                        "text/csv",
                        help="Download primer pairs as CSV table"
                    )
                with col2:
                    st.download_button(
                        "📥 Download Primer Sequences (FASTA)",
                        "\n".join(fasta_out),
                        "primers.fasta",
                        "text/plain",
                        help="Download primer sequences in FASTA format"
                    )
                
                # Visualization
                st.plotly_chart(
                    plot_primers(len(st.session_state.seq), st.session_state.pairs, st.session_state.params['primer_len']),
                    use_container_width=True
                )
        
        else:
            st.info("Please input a DNA sequence to design primers")
    
    # Footer
    st.markdown("---")
    st.markdown(
        "<div style='text-align: center; color: #666;'>"
        "SimplePrimer v1.0 | Automated PCR Primer Design"
        "</div>",
        unsafe_allow_html=True
    )
    
    # Reset button
    if st.button("🔄 Reset Application", help="Clear all data and start fresh"):
        reset_session()
        st.rerun()

if __name__ == "__main__":
    main()

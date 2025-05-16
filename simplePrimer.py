import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from io import StringIO
import pandas as pd

def is_dna(seq):
    return all(base in "ATCG" for base in seq)

def analyze_dna(seq):
    return gc_fraction(seq) * 100, mt.Tm_Wallace(seq)

def reverse_complement(seq):
    complement = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement)[::-1]

def design_primers(seq, primer_len=20):
    fwd = seq[:primer_len]
    rev = reverse_complement(seq[-primer_len:])
    return {
        "Forward Primer": fwd,
        "Forward GC%": gc_fraction(fwd) * 100,
        "Forward Tm": mt.Tm_Wallace(fwd),
        "Reverse Primer": rev,
        "Reverse GC%": gc_fraction(rev) * 100,
        "Reverse Tm": mt.Tm_Wallace(rev)
    }

st.title("Simple DNA Analyzer + Primer Designer")

# --- Text Input ---
st.subheader("Paste FASTA Sequence")
fasta_input = st.text_area("Paste your FASTA sequence here", height=150)

def parse_fasta(text):
    lines = text.strip().splitlines()
    if lines and lines[0].startswith(">"):
        return ''.join(lines[1:]).upper()
    return None

seq = None

# Check pasted input
if fasta_input:
    seq = parse_fasta(fasta_input)

# --- File Upload ---
st.subheader("Or Upload FASTA File")
uploaded_file = st.file_uploader("Upload .fasta, .fa, or .txt", type=["fasta", "fa", "txt"])
if uploaded_file and not seq:
    try:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        record = SeqIO.read(stringio, "fasta")
        seq = str(record.seq).upper()
    except Exception as e:
        st.error(f"File error: {e}")

# --- Analysis & Primer Output ---
if seq:
    if is_dna(seq):
        st.success("Valid DNA Sequence ✅")
        gc_content, tm = analyze_dna(seq)
        st.write(f"**GC Content:** {gc_content:.2f}%")
        st.write(f"**Melting Temp (Tm):** {tm:.2f}°C")

        primer_len = st.slider("Primer Length", 16, 25, 20)
        primers = design_primers(seq, primer_len)

        st.subheader("Designed Primers")
        st.write(f"**Forward Primer:** {primers['Forward Primer']}")
        st.write(f"GC%: {primers['Forward GC%']:.2f} | Tm: {primers['Forward Tm']:.2f}°C")

        st.write(f"**Reverse Primer:** {primers['Reverse Primer']}")
        st.write(f"GC%: {primers['Reverse GC%']:.2f} | Tm: {primers['Reverse Tm']:.2f}°C")

        df = pd.DataFrame([primers])
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Primers as CSV", csv, "primers.csv", "text/csv")
    else:
        st.error("Invalid DNA sequence ❌")
else:
    st.info("Paste or upload a FASTA sequence to start.")

import gradio as gr
import matplotlib.pyplot as plt
import pandas as pd
import io
from Bio import AlignIO, Phylo
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- 1. BIOPHYSICAL & EVOLUTIONARY LOGIC ---

def analyze_protein_biophysics(aln):
    """Calculates research-grade biochemical parameters of the reference sequence."""
    # Using the first sequence as the reference
    ref_seq = str(aln[0].seq).replace("-", "").upper()
    # Filter for standard amino acids to ensure ProtParam compatibility
    clean_seq = "".join([aa for aa in ref_seq if aa in "ACDEFGHIKLMNPQRSTVWY"])
    analysed = ProteinAnalysis(clean_seq)
    
    return {
        "Molecular Weight (Da)": round(analysed.molecular_weight(), 2),
        "Isoelectric Point (pI)": round(analysed.isoelectric_point(), 2),
        "Instability Index": round(analysed.instability_index(), 2),
        "Aromaticity": round(analysed.aromaticity(), 4)
    }

def calculate_conservation(aln):
    """Calculates column-wise identity scores across the alignment."""
    scores = []
    for i in range(aln.get_alignment_length()):
        column = aln[:, i]
        filtered = column.replace("-", "")
        if not filtered:
            scores.append(0); continue
        counts = [column.count(res) for res in set(filtered)]
        scores.append(max(counts) / len(aln))
    return scores

# --- 2. MAIN PROCESSING PIPELINE ---

def full_bio_pipeline(file):
    try:
        # Load Alignment (Detects Clustal or FASTA format)
        try:
            alignment = AlignIO.read(file.name, "clustal")
        except:
            alignment = AlignIO.read(file.name, "fasta")

        # 1. Biophysical Parameters
        bio_params = analyze_protein_biophysics(alignment)
        params_df = pd.DataFrame([bio_params])
        
        # 2. Alignment Metadata Stats
        meta_df = pd.DataFrame([{
            "Total Sequences": len(alignment),
            "Alignment Length": alignment.get_alignment_length(),
            "Gap Frequency (%)": round((sum(r.seq.count('-') for r in alignment)/(len(alignment)*len(alignment[0])))*100, 1)
        }])

        # 3. Evolutionary Conservation Map
        scores = calculate_conservation(alignment)
        fig_cons, ax_c = plt.subplots(figsize=(10, 4))
        ax_c.plot(scores, color='#1e3d59', linewidth=1)
        ax_c.fill_between(range(len(scores)), scores, color='#3498db', alpha=0.2)
        ax_c.axhline(0.9, color='red', linestyle='--', alpha=0.5, label="Functional Threshold (90%)")
        ax_c.set_title("Sequence Conservation & Evolutionary Constraint")
        ax_c.set_xlabel("Residue Position")
        ax_c.set_ylabel("Identity Score")
        ax_c.legend()
        plt.tight_layout()

        # 4. Phylogenetic Tree (NJ Method)
        calc = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calc, 'nj')
        tree = constructor.build_tree(alignment)
        
        fig_tree = plt.figure(figsize=(8, 6))
        ax_t = fig_tree.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax_t, do_show=False)
        ax_t.set_title("Evolutionary Lineage (Neighbor-Joining)")
        plt.tight_layout()

        return fig_cons, fig_tree, params_df, meta_df, "Analysis Complete!"

    except Exception as e:
        return None, None, None, None, f"Error: {str(e)}"

# --- 3. GRADIO INTERFACE ---

with gr.Blocks(theme=gr.themes.Soft(primary_hue="blue")) as demo:
    gr.Markdown("# 🧬 Bio-Insight PRO: End-to-End Molecular Suite")
    gr.Markdown("An integrated suite for sequence identity, phylogeny, and biophysical analysis.")
    
    with gr.Row():
        with gr.Column(scale=1):
            file_input = gr.File(label="Upload Alignment (.aln or .fasta)")
            run_btn = gr.Button("🚀 Execute Deep Analysis", variant="primary")
            status = gr.Textbox(label="Status", interactive=False)
            
            gr.Markdown("### 📊 Biophysical Properties")
            biophys_table = gr.DataFrame()
            
            gr.Markdown("### ⚙️ Alignment Metadata")
            meta_table = gr.DataFrame()

        with gr.Column(scale=2):
            with gr.Tabs():
                with gr.TabItem("Evolutionary Map"):
                    plot_cons = gr.Plot()
                    gr.Markdown("**Interpretation:** Higher peaks indicate regions of high evolutionary constraint, representing critical functional domains.")
                
                with gr.TabItem("Phylogenetic Tree"):
                    plot_tree = gr.Plot()
                    gr.Markdown("**Note:** Branch lengths represent the estimated evolutionary distance between species.")

    run_btn.click(
        fn=full_bio_pipeline, 
        inputs=[file_input], 
        outputs=[plot_cons, plot_tree, biophys_table, meta_table, status]
    )

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 10000))
    # show_api=False is the key to stopping the TypeError
    demo.launch(
        server_name="0.0.0.0", 
        server_port=port, 
        show_api=False
    )

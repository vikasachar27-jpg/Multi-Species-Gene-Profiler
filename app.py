import gradio as gr
import matplotlib.pyplot as plt
import io
import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- UNIQUE ANALYTICAL FUNCTIONS ---

def get_alignment_stats(aln):
    """Generates unique metadata about the alignment."""
    stats = {
        "Total Sequences": len(aln),
        "Alignment Length": aln.get_alignment_length(),
        "Avg. Gap Density (%)": round((sum(rec.seq.count("-") for rec in aln) / (len(aln) * aln.get_alignment_length())) * 100, 2)
    }
    return pd.DataFrame([stats])

def calculate_gc_content(aln):
    """Unique Feature: Calculates GC content across the aligned sequences."""
    gc_scores = []
    for record in aln:
        s = str(record.seq).upper()
        gc = (s.count('G') + s.count('C')) / len(s) * 100 if len(s) > 0 else 0
        gc_scores.append(round(gc, 2))
    return gc_scores

def calculate_conservation(aln):
    """Calculates residue identity scores."""
    scores = []
    length = aln.get_alignment_length()
    for i in range(length):
        column = aln[:, i]
        filtered = column.replace("-", "")
        if not filtered:
            scores.append(0); continue
        counts = [column.count(res) for res in set(filtered)]
        scores.append(max(counts) / len(aln))
    return scores

# --- MAIN PIPELINE ---

def process_and_analyze(file):
    try:
        # 1. Load Alignment
        try:
            alignment = AlignIO.read(file.name, "clustal")
        except:
            alignment = AlignIO.read(file.name, "fasta")

        # 2. Generate Basic Graphs
        scores = calculate_conservation(alignment)
        fig_cons, ax = plt.subplots(figsize=(10, 4))
        ax.plot(scores, color='#2c3e50', linewidth=1)
        ax.fill_between(range(len(scores)), scores, color='#3498db', alpha=0.3)
        ax.set_title("Functional Conservation Profile", fontsize=12)
        plt.tight_layout()

        # 3. Generate Phylogenetic Tree
        calc = DistanceCalculator('identity')
        tree = DistanceTreeConstructor(calc, 'nj').build_tree(alignment)
        fig_tree = plt.figure(figsize=(8, 6))
        ax_t = fig_tree.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax_t, do_show=False)
        plt.tight_layout()

        # 4. UNIQUE DATA: GC Content & Stats
        gc_data = calculate_gc_content(alignment)
        species_names = [r.id for r in alignment]
        fig_gc, ax_gc = plt.subplots(figsize=(8, 4))
        ax_gc.bar(species_names, gc_data, color=['#e67e22', '#27ae60', '#8e44ad', '#c0392b'])
        ax_gc.set_title("Genomic GC Content Signature (%)")
        plt.xticks(rotation=45)
        plt.tight_layout()

        stats_df = get_alignment_stats(alignment)

        return fig_cons, fig_tree, fig_gc, stats_df, "Success: Bio-Analysis Complete!"

    except Exception as e:
        return None, None, None, None, f"Error: {str(e)}"

# --- GRADIO INTERFACE ---

with gr.Blocks(theme=gr.themes.Soft(primary_hue="teal")) as demo:
    gr.Markdown("# 🧬 Bio-Evolutionary Advanced Dashboard")
    gr.Markdown("Upload an MSA file to explore conservation, phylogeny, and unique genomic signatures.")
    
    with gr.Row():
        with gr.Column(scale=1):
            file_input = gr.File(label="Upload .aln or .fasta Alignment")
            run_btn = gr.Button("🔍 Deep Bio-Analysis", variant="primary")
            status = gr.Textbox(label="Status",lines=6)
            stats_table = gr.DataFrame(label="Alignment Metadata")

        with gr.Column(scale=2):
            with gr.Tabs():
                with gr.TabItem("Conservation Map"):
                    plot_cons = gr.Plot()
                with gr.TabItem("Phylogenetic Tree"):
                    plot_tree = gr.Plot()
                with gr.TabItem("GC Signature"):
                    plot_gc = gr.Plot()

    run_btn.click(
        fn=process_and_analyze, 
        inputs=file_input, 
        outputs=[plot_cons, plot_tree, plot_gc, stats_table, status]
    )

if __name__ == "__main__":
    demo.launch()
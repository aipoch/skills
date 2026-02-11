#!/usr/bin/env python3
"""
GO/KEGG Enrichment Analysis Pipeline (Pure Python)

Automated pipeline for Gene Ontology and KEGG pathway enrichment analysis.
Uses gseapy (Python) instead of clusterProfiler (R).
Supports multiple organisms, ID types, and generates comprehensive visualizations.

Author: AI Assistant
Technical Difficulty: High
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np

# Try to import gseapy
try:
    import gseapy as gp
    from gseapy.plot import barplot, dotplot
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    print("⚠️  Warning: gseapy not installed. Install with: pip install gseapy")

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Organism mapping configuration
ORGANISM_CONFIG = {
    "human": {
        "scientific_name": "Homo sapiens",
        "kegg_code": "hsa",
        "tax_id": "9606",
        "enrichr_library": "Human"
    },
    "mouse": {
        "scientific_name": "Mus musculus", 
        "kegg_code": "mmu",
        "tax_id": "10090",
        "enrichr_library": "Mouse"
    },
    "rat": {
        "scientific_name": "Rattus norvegicus",
        "kegg_code": "rno", 
        "tax_id": "10116",
        "enrichr_library": "Rat"
    },
    "zebrafish": {
        "scientific_name": "Danio rerio",
        "kegg_code": "dre",
        "tax_id": "7955",
        "enrichr_library": "Zebrafish"
    },
    "fly": {
        "scientific_name": "Drosophila melanogaster",
        "kegg_code": "dme",
        "tax_id": "7227",
        "enrichr_library": "Fly"
    },
    "yeast": {
        "scientific_name": "Saccharomyces cerevisiae",
        "kegg_code": "sce",
        "tax_id": "4932",
        "enrichr_library": "Yeast"
    }
}

# GO Enrichr library mapping
GO_LIBRARIES = {
    "BP": "GO_Biological_Process_2021",
    "MF": "GO_Molecular_Function_2021",
    "CC": "GO_Cellular_Component_2021"
}


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="GO/KEGG Enrichment Analysis Pipeline (Pure Python)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic GO and KEGG analysis
  python main.py --genes gene_list.txt --organism human
  
  # KEGG only with custom cutoff
  python main.py --genes genes.txt --analysis kegg --pvalue-cutoff 0.01
  
  # GO only with specific ontologies
  python main.py --genes genes.txt --analysis go --go-ontologies BP,MF
  
  # Use Enrichr (online, faster)
  python main.py --genes genes.txt --use-enrichr --organism human
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--genes", "-g",
        type=str,
        required=True,
        help="Path to gene list file (one gene per line) or comma-separated gene list"
    )
    
    # Optional arguments
    parser.add_argument(
        "--organism", "-o",
        type=str,
        default="human",
        choices=list(ORGANISM_CONFIG.keys()),
        help="Organism (default: human)"
    )
    parser.add_argument(
        "--id-type",
        type=str,
        default="symbol",
        choices=["symbol", "entrez", "ensembl", "refseq"],
        help="Gene ID type (default: symbol)"
    )
    parser.add_argument(
        "--background", "-b",
        type=str,
        default=None,
        help="Background gene list file (default: all genes)"
    )
    parser.add_argument(
        "--analysis", "-a",
        type=str,
        default="all",
        choices=["go", "kegg", "all"],
        help="Analysis type (default: all)"
    )
    parser.add_argument(
        "--go-ontologies",
        type=str,
        default="BP,MF,CC",
        help="GO ontologies to analyze, comma-separated (default: BP,MF,CC)"
    )
    parser.add_argument(
        "--pvalue-cutoff",
        type=float,
        default=0.05,
        help="P-value cutoff for significance (default: 0.05)"
    )
    parser.add_argument(
        "--qvalue-cutoff",
        type=float,
        default=0.2,
        help="Adjusted p-value (q-value) cutoff (default: 0.2)"
    )
    parser.add_argument(
        "--output", "-out",
        type=str,
        default="./enrichment_results",
        help="Output directory (default: ./enrichment_results)"
    )
    parser.add_argument(
        "--format",
        type=str,
        default="csv",
        choices=["csv", "tsv", "excel", "all"],
        help="Output format (default: csv)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top terms to visualize (default: 20)"
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=5,
        help="Minimum genes in category (default: 5)"
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=500,
        help="Maximum genes in category (default: 500)"
    )
    parser.add_argument(
        "--use-enrichr",
        action="store_true",
        help="Use Enrichr API instead of local gseapy (faster, no download needed)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()


def validate_inputs(args: argparse.Namespace) -> Tuple[bool, str]:
    """Validate input files and parameters."""
    # Check if gseapy is available
    if not GSEAPY_AVAILABLE and not args.use_enrichr:
        return False, "gseapy not installed. Install with: pip install gseapy"
    
    # Check gene file exists or is a comma-separated list
    if not os.path.exists(args.genes) and ',' not in args.genes:
        return False, f"Gene list file not found: {args.genes}"
    
    # Check background file if provided
    if args.background and not os.path.exists(args.background):
        return False, f"Background file not found: {args.background}"
    
    # Check cutoffs are valid
    if not 0 < args.pvalue_cutoff <= 1:
        return False, "P-value cutoff must be between 0 and 1"
    if not 0 < args.qvalue_cutoff <= 1:
        return False, "Q-value cutoff must be between 0 and 1"
    
    return True, "Validation passed"


def read_gene_list(input_data: str) -> List[str]:
    """Read gene list from file or string."""
    genes = []
    
    # Check if it's a file path or comma-separated list
    if os.path.exists(input_data):
        # Read from file
        with open(input_data, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Handle CSV/TSV format
                    if ',' in line:
                        genes.append(line.split(',')[0])
                    elif '\t' in line:
                        genes.append(line.split('\t')[0])
                    else:
                        genes.append(line)
    elif ',' in input_data:
        # Comma-separated list
        genes = [g.strip() for g in input_data.split(',') if g.strip()]
    else:
        # Single gene
        genes = [input_data.strip()]
    
    return list(set(genes))  # Remove duplicates


def create_output_directories(base_path: str) -> Dict[str, str]:
    """Create output directory structure."""
    dirs = {
        'base': base_path,
        'go': os.path.join(base_path, 'go_enrichment'),
        'kegg': os.path.join(base_path, 'kegg_enrichment'),
        'visualization': os.path.join(base_path, 'visualization')
    }
    
    for path in dirs.values():
        os.makedirs(path, exist_ok=True)
    
    return dirs


def run_go_enrichment_local(
    gene_list: List[str],
    organism: str,
    ontology: str,
    output_dir: str,
    args: argparse.Namespace
) -> Optional[pd.DataFrame]:
    """Run GO enrichment analysis using gseapy locally."""
    
    print(f"\n  Processing {ontology}...")
    
    try:
        # Map organism to gseapy organism
        organism_map = {
            "human": "human",
            "mouse": "mouse",
            "rat": "rat",
            "zebrafish": "zebrafish",
            "fly": "fly",
            "yeast": "yeast"
        }
        
        gsea_organism = organism_map.get(organism, "human")
        
        # Run enrichment
        enr = gp.enrichGO(
            gene_list=gene_list,
            organism=gsea_organism,
            geneid=args.id_type.upper() if args.id_type != 'symbol' else 'SYMBOL',
            ont=ontology,
            pvalue_cutoff=args.pvalue_cutoff,
            qvalue_cutoff=args.qvalue_cutoff,
            min_gene_size=args.min_genes,
            max_gene_size=args.max_genes,
            background=args.background if args.background else None,
            outdir=None,  # We'll handle output ourselves
            verbose=args.verbose
        )
        
        if enr is None or enr.res2d is None or len(enr.res2d) == 0:
            print(f"    No significant {ontology} terms found")
            return None
        
        # Get results
        results = enr.res2d.copy()
        
        # Save results
        output_file = os.path.join(output_dir, f"GO_{ontology}_results.csv")
        
        # Format output
        if args.format in ["csv", "all"]:
            results.to_csv(output_file, index=False)
            print(f"    Saved: {output_file}")
        
        if args.format in ["tsv", "all"]:
            tsv_file = os.path.join(output_dir, f"GO_{ontology}_results.tsv")
            results.to_csv(tsv_file, sep='\t', index=False)
        
        if args.format in ["excel", "all"]:
            excel_file = os.path.join(output_dir, f"GO_{ontology}_results.xlsx")
            results.to_excel(excel_file, index=False)
        
        # Create visualizations
        if MATPLOTLIB_AVAILABLE and len(results) > 0:
            create_visualizations(results, ontology, "GO", output_dir, args)
        
        return results
        
    except Exception as e:
        print(f"    Error in GO {ontology}: {str(e)}")
        return None


def run_kegg_enrichment_local(
    gene_list: List[str],
    organism: str,
    output_dir: str,
    args: argparse.Namespace
) -> Optional[pd.DataFrame]:
    """Run KEGG enrichment analysis using gseapy locally."""
    
    print(f"\n  Processing KEGG pathways...")
    
    try:
        # Get organism code
        kegg_code = ORGANISM_CONFIG[organism]["kegg_code"]
        
        # Run enrichment
        enr = gp.enrichKEGG(
            gene_list=gene_list,
            organism=kegg_code,
            geneid=args.id_type.upper() if args.id_type != 'symbol' else 'SYMBOL',
            pvalue_cutoff=args.pvalue_cutoff,
            qvalue_cutoff=args.qvalue_cutoff,
            min_gene_size=args.min_genes,
            max_gene_size=args.max_genes,
            background=args.background if args.background else None,
            outdir=None,
            verbose=args.verbose
        )
        
        if enr is None or enr.res2d is None or len(enr.res2d) == 0:
            print(f"    No significant KEGG pathways found")
            return None
        
        # Get results
        results = enr.res2d.copy()
        
        # Save results
        output_file = os.path.join(output_dir, "KEGG_results.csv")
        
        # Format output
        if args.format in ["csv", "all"]:
            results.to_csv(output_file, index=False)
            print(f"    Saved: {output_file}")
        
        if args.format in ["tsv", "all"]:
            tsv_file = os.path.join(output_dir, "KEGG_results.tsv")
            results.to_csv(tsv_file, sep='\t', index=False)
        
        if args.format in ["excel", "all"]:
            excel_file = os.path.join(output_dir, "KEGG_results.xlsx")
            results.to_excel(excel_file, index=False)
        
        # Create visualizations
        if MATPLOTLIB_AVAILABLE and len(results) > 0:
            create_visualizations(results, "Pathway", "KEGG", output_dir, args)
        
        return results
        
    except Exception as e:
        print(f"    Error in KEGG: {str(e)}")
        return None


def run_enrichr_analysis(
    gene_list: List[str],
    organism: str,
    analysis_type: str,
    go_ontologies: List[str],
    output_dirs: Dict[str, str],
    args: argparse.Namespace
) -> Dict[str, pd.DataFrame]:
    """Run enrichment analysis using Enrichr API."""
    
    results_summary = {}
    
    # Map organism
    organism_map = {
        "human": "Human",
        "mouse": "Mouse",
        "rat": "Rat",
        "zebrafish": "Zebrafish",
        "fly": "Fly",
        "yeast": "Yeast"
    }
    
    enrichr_organism = organism_map.get(organism, "Human")
    
    # GO Enrichment
    if analysis_type in ["go", "all"]:
        print("\n=== Running GO Enrichment Analysis (Enrichr) ===")
        
        for ontology in go_ontologies:
            library = GO_LIBRARIES.get(ontology.upper())
            if not library:
                continue
            
            print(f"\n  Processing {ontology}...")
            
            try:
                # Run Enrichr
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=library,
                    organism=enrichr_organism,
                    outdir=None,
                    cutoff=args.pvalue_cutoff
                )
                
                if enr is None or enr.res2d is None or len(enr.res2d) == 0:
                    print(f"    No significant {ontology} terms found")
                    results_summary[f"GO_{ontology}"] = None
                    continue
                
                # Get results
                results = enr.res2d.copy()
                
                # Filter by q-value
                if 'Adjusted P-value' in results.columns:
                    results = results[results['Adjusted P-value'] <= args.qvalue_cutoff]
                
                if len(results) == 0:
                    print(f"    No significant {ontology} terms after filtering")
                    results_summary[f"GO_{ontology}"] = None
                    continue
                
                # Save results
                output_file = os.path.join(output_dirs['go'], f"GO_{ontology}_results.csv")
                
                if args.format in ["csv", "all"]:
                    results.to_csv(output_file, index=False)
                    print(f"    Saved: {output_file}")
                
                if args.format in ["tsv", "all"]:
                    results.to_csv(output_file.replace('.csv', '.tsv'), sep='\t', index=False)
                
                # Create visualizations
                if MATPLOTLIB_AVAILABLE:
                    create_enrichr_visualizations(results, ontology, output_dirs['go'], args)
                
                results_summary[f"GO_{ontology}"] = results
                
            except Exception as e:
                print(f"    Error in GO {ontology}: {str(e)}")
                results_summary[f"GO_{ontology}"] = None
    
    # KEGG Enrichment
    if analysis_type in ["kegg", "all"]:
        print("\n=== Running KEGG Enrichment Analysis (Enrichr) ===")
        
        try:
            # KEGG library names vary by organism
            kegg_library = f"KEGG_{enrichr_organism}_2019"
            
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=kegg_library,
                organism=enrichr_organism,
                outdir=None,
                cutoff=args.pvalue_cutoff
            )
            
            if enr is None or enr.res2d is None or len(enr.res2d) == 0:
                print(f"    No significant KEGG pathways found")
                results_summary["KEGG"] = None
            else:
                results = enr.res2d.copy()
                
                # Filter by q-value
                if 'Adjusted P-value' in results.columns:
                    results = results[results['Adjusted P-value'] <= args.qvalue_cutoff]
                
                if len(results) == 0:
                    print(f"    No significant KEGG pathways after filtering")
                    results_summary["KEGG"] = None
                else:
                    # Save results
                    output_file = os.path.join(output_dirs['kegg'], "KEGG_results.csv")
                    
                    if args.format in ["csv", "all"]:
                        results.to_csv(output_file, index=False)
                        print(f"    Saved: {output_file}")
                    
                    if args.format in ["tsv", "all"]:
                        results.to_csv(output_file.replace('.csv', '.tsv'), sep='\t', index=False)
                    
                    # Create visualizations
                    if MATPLOTLIB_AVAILABLE:
                        create_enrichr_visualizations(results, "Pathway", output_dirs['kegg'], args)
                    
                    results_summary["KEGG"] = results
                    
        except Exception as e:
            print(f"    Error in KEGG: {str(e)}")
            results_summary["KEGG"] = None
    
    return results_summary


def create_visualizations(
    results: pd.DataFrame,
    category: str,
    analysis_type: str,
    output_dir: str,
    args: argparse.Namespace
):
    """Create visualization plots."""
    
    if not MATPLOTLIB_AVAILABLE:
        return
    
    try:
        # Prepare data for plotting
        plot_data = results.head(args.top_n).copy()
        
        if len(plot_data) == 0:
            return
        
        # Bar plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Use Term or Description column
        term_col = 'Term' if 'Term' in plot_data.columns else 'Description'
        pval_col = 'pvalue' if 'pvalue' in plot_data.columns else 'P-value'
        
        y_pos = np.arange(len(plot_data))
        
        # Color by p-value
        colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(plot_data)))
        
        ax.barh(y_pos, -np.log10(plot_data[pval_col]), color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data[term_col], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel('-log10(p-value)', fontsize=10)
        ax.set_title(f'{analysis_type} {category} Enrichment', fontsize=12)
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        barplot_file = os.path.join(output_dir, f"{analysis_type}_{category}_barplot.png")
        plt.savefig(barplot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    Created barplot: {barplot_file}")
        
        # Dot plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get gene ratio
        if 'GeneRatio' in plot_data.columns:
            ratio_col = 'GeneRatio'
        elif 'Overlap' in plot_data.columns:
            ratio_col = 'Overlap'
        else:
            ratio_col = None
        
        if ratio_col:
            # Convert ratio to numeric if needed
            if plot_data[ratio_col].dtype == object:
                plot_data['ratio_numeric'] = plot_data[ratio_col].apply(
                    lambda x: eval(x) if '/' in str(x) else float(x)
                )
            else:
                plot_data['ratio_numeric'] = plot_data[ratio_col]
            
            scatter = ax.scatter(
                plot_data['ratio_numeric'],
                y_pos,
                s=100,
                c=-np.log10(plot_data[pval_col]),
                cmap='RdYlBu_r',
                alpha=0.7
            )
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels(plot_data[term_col], fontsize=8)
            ax.invert_yaxis()
            ax.set_xlabel('Gene Ratio', fontsize=10)
            ax.set_title(f'{analysis_type} {category} Enrichment (Dot Plot)', fontsize=12)
            
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('-log10(p-value)', rotation=270, labelpad=15)
            
            plt.tight_layout()
            dotplot_file = os.path.join(output_dir, f"{analysis_type}_{category}_dotplot.png")
            plt.savefig(dotplot_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    Created dotplot: {dotplot_file}")
        
    except Exception as e:
        print(f"    Warning: Could not create visualizations: {e}")


def create_enrichr_visualizations(
    results: pd.DataFrame,
    category: str,
    output_dir: str,
    args: argparse.Namespace
):
    """Create visualizations for Enrichr results."""
    
    if not MATPLOTLIB_AVAILABLE:
        return
    
    try:
        # Prepare data
        plot_data = results.head(args.top_n).copy()
        
        if len(plot_data) == 0:
            return
        
        # Get column names
        term_col = 'Term' if 'Term' in plot_data.columns else 'Description'
        pval_col = 'P-value' if 'P-value' in plot_data.columns else 'Adjusted P-value'
        
        # Bar plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        y_pos = np.arange(len(plot_data))
        colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(plot_data)))
        
        ax.barh(y_pos, -np.log10(plot_data[pval_col]), color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data[term_col].str[:50], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel('-log10(p-value)', fontsize=10)
        ax.set_title(f'{category} Enrichment', fontsize=12)
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        barplot_file = os.path.join(output_dir, f"{category}_barplot.png")
        plt.savefig(barplot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    Created barplot: {barplot_file}")
        
    except Exception as e:
        print(f"    Warning: Could not create visualizations: {e}")


def generate_summary_report(
    results_summary: Dict[str, Optional[pd.DataFrame]],
    output_dirs: Dict[str, str],
    args: argparse.Namespace
) -> str:
    """Generate human-readable summary report."""
    
    report_lines = [
        "=" * 60,
        "GO/KEGG ENRICHMENT ANALYSIS REPORT",
        "=" * 60,
        "",
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Analysis Method: {'Enrichr API' if args.use_enrichr else 'Local gseapy'}",
        f"Organism: {ORGANISM_CONFIG[args.organism]['scientific_name']}",
        f"Gene ID Type: {args.id_type}",
        f"Analysis Type: {args.analysis}",
        f"P-value Cutoff: {args.pvalue_cutoff}",
        f"Q-value Cutoff: {args.qvalue_cutoff}",
        "",
        "-" * 60,
        "ENRICHMENT RESULTS SUMMARY",
        "-" * 60,
        "",
    ]
    
    # Summarize results
    for analysis_name, results in results_summary.items():
        if results is not None and len(results) > 0:
            report_lines.append(f"{analysis_name}:")
            report_lines.append(f"  - Significant terms: {len(results)}")
            
            # Show top terms
            term_col = 'Term' if 'Term' in results.columns else 'Description'
            top_terms = results[term_col].head(5).tolist()
            report_lines.append(f"  - Top terms: {', '.join(top_terms)}")
            report_lines.append("")
        else:
            report_lines.append(f"{analysis_name}: No significant results")
            report_lines.append("")
    
    report_lines.extend([
        "-" * 60,
        "OUTPUT FILES",
        "-" * 60,
        "",
    ])
    
    # List generated files
    if os.path.exists(output_dirs['base']):
        for root, dirs, files in os.walk(output_dirs['base']):
            level = root.replace(output_dirs['base'], '').count(os.sep)
            indent = '  ' * level
            subdir = os.path.basename(root)
            if subdir and subdir != os.path.basename(output_dirs['base']):
                report_lines.append(f"{indent}{subdir}/")
            
            subindent = '  ' * (level + 1)
            for file in sorted(files):
                if not file.startswith('.'):
                    report_lines.append(f"{subindent}{file}")
    
    report_lines.extend([
        "",
        "-" * 60,
        "INTERPRETATION GUIDE",
        "-" * 60,
        "",
        "1. GO Enrichment Results:",
        "   - BP: Biological Process - what biological processes are affected",
        "   - MF: Molecular Function - what molecular functions are involved", 
        "   - CC: Cellular Component - where in the cell the genes act",
        "",
        "2. KEGG Pathway Results:",
        "   - Shows affected metabolic and signaling pathways",
        "",
        "3. Key Statistics:",
        "   - GeneRatio: Proportion of input genes in the term",
        "   - P-value: Statistical significance",
        "   - Adjusted P-value: Benjamini-Hochberg corrected p-value",
        "",
        "=" * 60,
    ])
    
    return '\n'.join(report_lines)


def main():
    """Main entry point."""
    print("=" * 60)
    print("GO/KEGG Enrichment Analysis Pipeline (Pure Python)")
    print("=" * 60)
    
    # Parse arguments
    args = parse_arguments()
    
    # Validate inputs
    valid, message = validate_inputs(args)
    if not valid:
        print(f"ERROR: {message}")
        sys.exit(1)
    
    if args.verbose:
        print(f"Validation: {message}")
    
    # Read gene list
    print(f"\nReading gene list from: {args.genes}")
    genes = read_gene_list(args.genes)
    print(f"Loaded {len(genes)} unique genes")
    
    if len(genes) == 0:
        print("ERROR: No genes found in input file")
        sys.exit(1)
    
    if len(genes) < 5:
        print("WARNING: Very few genes provided. Results may be limited.")
    
    if args.verbose:
        print(f"First 10 genes: {', '.join(genes[:10])}")
    
    # Create output directories
    print(f"\nCreating output directory: {args.output}")
    output_dirs = create_output_directories(args.output)
    
    # Parse GO ontologies
    go_ontologies = [o.strip().upper() for o in args.go_ontologies.split(',')]
    
    # Run enrichment analysis
    results_summary = {}
    
    if args.use_enrichr:
        # Use Enrichr API
        print("\nRunning enrichment analysis using Enrichr API...")
        results_summary = run_enrichr_analysis(
            genes, args.organism, args.analysis, go_ontologies, output_dirs, args
        )
    else:
        # Use local gseapy
        print("\nRunning enrichment analysis using local gseapy...")
        
        # GO Enrichment
        if args.analysis in ["go", "all"]:
            print("\n=== Running GO Enrichment Analysis ===")
            for ontology in go_ontologies:
                results = run_go_enrichment_local(
                    genes, args.organism, ontology, output_dirs['go'], args
                )
                results_summary[f"GO_{ontology}"] = results
        
        # KEGG Enrichment
        if args.analysis in ["kegg", "all"]:
            print("\n=== Running KEGG Enrichment Analysis ===")
            results = run_kegg_enrichment_local(
                genes, args.organism, output_dirs['kegg'], args
            )
            results_summary["KEGG"] = results
    
    # Generate summary report
    print("\nGenerating summary report...")
    report = generate_summary_report(results_summary, output_dirs, args)
    report_path = os.path.join(output_dirs['base'], 'REPORT.txt')
    with open(report_path, 'w') as f:
        f.write(report)
    
    print(f"\nReport saved to: {report_path}")
    print("\n" + report)
    
    # Success message
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nResults available in: {os.path.abspath(args.output)}")
    print("\nNext steps:")
    print("  1. Review CSV files for detailed statistics")
    print("  2. Check PNG visualizations in output folders")
    print("  3. Read REPORT.txt for interpretation guidance")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

"""
generate_consolidated_report.py - Generate MultiQC-style consolidated reports
Usage: python generate_consolidated_report.py <results_dir> <output_prefix> [--format html,txt]
"""

import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import json
import base64
from io import BytesIO

def load_all_results(results_dir):
    """Load all results from the pipeline output."""
    results_dir = Path(results_dir)
    
    data = {}
    
    # Load copy number results
    cn_file = results_dir / 'cnv_calls' / 'copy_numbers_gene_level.txt'
    if cn_file.exists():
        data['gene_cn'] = pd.read_csv(cn_file, sep='\t')
    
    cn_exon_file = results_dir / 'cnv_calls' / 'copy_numbers.txt'
    if cn_exon_file.exists():
        data['exon_cn'] = pd.read_csv(cn_exon_file, sep='\t')
    
    # Load Z-scores
    z_scores_file = results_dir / 'normalized' / 'z_scores.txt'
    if z_scores_file.exists():
        data['z_scores'] = pd.read_csv(z_scores_file, sep='\t')
    
    # Load reference statistics
    ref_stats_file = results_dir / 'normalized' / 'z_scores_ref_stats.txt'
    if ref_stats_file.exists():
        data['ref_stats'] = pd.read_csv(ref_stats_file, sep='\t')
    
    # Load allele counts
    allele_file = results_dir / 'allele_counts' / 'allele_counts.txt'
    if allele_file.exists():
        data['alleles'] = pd.read_csv(allele_file, sep='\t')
    
    # Load sample info
    sample_info_file = results_dir / 'allele_counts' / 'sample_info.txt'
    if sample_info_file.exists():
        data['sample_info'] = pd.read_csv(sample_info_file, sep='\t')
    
    # Load bootstrap results if available
    bootstrap_file = results_dir / 'cnv_calls' / 'copy_numbers_bootstrap.txt'
    if bootstrap_file.exists():
        data['bootstrap'] = pd.read_csv(bootstrap_file, sep='\t')
    
    return data

def generate_summary_statistics(data):
    """Generate comprehensive summary statistics."""
    summary = {
        'run_info': {
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_samples': 0,
            'reference_samples': 0,
            'test_samples': 0,
            'exons_analyzed': [],
            'genes_analyzed': []
        },
        'quality_metrics': {},
        'cnv_summary': {},
        'clinical_summary': {}
    }
    
    # Basic run information
    if 'sample_info' in data:
        sample_info = data['sample_info']
        summary['run_info']['total_samples'] = len(sample_info)
        summary['run_info']['reference_samples'] = len(sample_info[sample_info['sample_type'] == 'reference'])
        summary['run_info']['test_samples'] = len(sample_info[sample_info['sample_type'] == 'test'])
    
    if 'exon_cn' in data:
        summary['run_info']['exons_analyzed'] = sorted(data['exon_cn']['exon'].unique().tolist())
    
    if 'gene_cn' in data:
        summary['run_info']['genes_analyzed'] = sorted(data['gene_cn']['gene'].unique().tolist())
    
    # Quality metrics
    if 'z_scores' in data:
        z_scores = data['z_scores']
        summary['quality_metrics'] = {
            'mean_z_score': float(z_scores['z_score'].mean()),
            'std_z_score': float(z_scores['z_score'].std()),
            'high_confidence_calls': int((z_scores.get('confidence_score', pd.Series([0.5])) > 0.8).sum()) if 'confidence_score' in z_scores.columns else 0,
            'low_outlier_probability_calls': int((z_scores.get('outlier_probability', pd.Series([0.5])) < 0.1).sum()) if 'outlier_probability' in z_scores.columns else 0
        }
    
    # CNV summary
    if 'gene_cn' in data:
        gene_cn = data['gene_cn']
        
        # SMN1 specific analysis
        smn1_data = gene_cn[gene_cn['gene'] == 'SMN1']
        if not smn1_data.empty:
            cn_dist = smn1_data['estimated_copy_number'].value_counts().sort_index()
            summary['cnv_summary']['SMN1'] = {
                'copy_number_distribution': cn_dist.to_dict(),
                'mean_copy_number': float(smn1_data['estimated_copy_number'].mean()),
                'affected_samples': int((smn1_data['estimated_copy_number'] <= 0.5).sum()),
                'carrier_samples': int(((smn1_data['estimated_copy_number'] > 0.5) & 
                                      (smn1_data['estimated_copy_number'] <= 1.5)).sum()),
                'normal_samples': int(((smn1_data['estimated_copy_number'] > 1.5) & 
                                     (smn1_data['estimated_copy_number'] <= 2.5)).sum()),
                'duplication_samples': int((smn1_data['estimated_copy_number'] > 2.5).sum())
            }
        
        # SMN2 specific analysis
        smn2_data = gene_cn[gene_cn['gene'] == 'SMN2']
        if not smn2_data.empty:
            cn_dist = smn2_data['estimated_copy_number'].value_counts().sort_index()
            summary['cnv_summary']['SMN2'] = {
                'copy_number_distribution': cn_dist.to_dict(),
                'mean_copy_number': float(smn2_data['estimated_copy_number'].mean())
            }
    
    # Clinical summary
    if 'gene_cn' in data:
        gene_cn = data['gene_cn']
        if 'clinical_significance' in gene_cn.columns:
            clin_sig = gene_cn[gene_cn['gene'] == 'SMN1']['clinical_significance'].value_counts()
            summary['clinical_summary'] = clin_sig.to_dict()
    
    return summary

def create_consolidated_plots(data, output_dir):
    """Create comprehensive plots for consolidated report."""
    plots = {}
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Copy number distribution overview
    if 'gene_cn' in data:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # SMN1 copy number distribution
        smn1_data = data['gene_cn'][data['gene_cn']['gene'] == 'SMN1']
        if not smn1_data.empty:
            cn_counts = smn1_data['estimated_copy_number'].round().value_counts().sort_index()
            colors = ['red' if cn <= 0.5 else 'orange' if cn <= 1.5 else 'green' if cn <= 2.5 else 'blue' 
                     for cn in cn_counts.index]
            
            axes[0,0].bar(cn_counts.index, cn_counts.values, color=colors, alpha=0.7)
            axes[0,0].set_title('SMN1 Copy Number Distribution')
            axes[0,0].set_xlabel('Copy Number')
            axes[0,0].set_ylabel('Sample Count')
            
            # Add clinical annotations
            axes[0,0].axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='SMA affected')
            axes[0,0].axvline(x=1.5, color='orange', linestyle='--', alpha=0.7, label='SMA carrier')
            axes[0,0].legend()
        
        # SMN2 copy number distribution
        smn2_data = data['gene_cn'][data['gene_cn']['gene'] == 'SMN2']
        if not smn2_data.empty:
            cn_counts = smn2_data['estimated_copy_number'].round().value_counts().sort_index()
            axes[0,1].bar(cn_counts.index, cn_counts.values, alpha=0.7, color='lightblue')
            axes[0,1].set_title('SMN2 Copy Number Distribution')
            axes[0,1].set_xlabel('Copy Number')
            axes[0,1].set_ylabel('Sample Count')
        
        # Confidence distribution
        if 'confidence_score' in data['gene_cn'].columns:
            conf_data = data['gene_cn']['confidence_score']
            axes[1,0].hist(conf_data, bins=20, alpha=0.7, edgecolor='black')
            axes[1,0].set_title('Confidence Score Distribution')
            axes[1,0].set_xlabel('Confidence Score')
            axes[1,0].set_ylabel('Sample Count')
            axes[1,0].axvline(x=0.8, color='red', linestyle='--', alpha=0.7, label='High confidence threshold')
            axes[1,0].legend()
        
        # Clinical significance pie chart
        if 'clinical_significance' in data['gene_cn'].columns:
            smn1_clin = data['gene_cn'][data['gene_cn']['gene'] == 'SMN1']['clinical_significance']
            clin_counts = smn1_clin.value_counts()
            
            colors = ['red', 'orange', 'green', 'blue', 'gray'][:len(clin_counts)]
            axes[1,1].pie(clin_counts.values, labels=clin_counts.index, autopct='%1.1f%%', 
                         colors=colors, startangle=90)
            axes[1,1].set_title('SMN1 Clinical Significance')
        
        plt.tight_layout()
        plot_file = plot_dir / 'overview_dashboard.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['overview_dashboard'] = plot_file
    
    # Plot 2: Quality metrics heatmap
    if 'z_scores' in data:
        z_scores = data['z_scores']
        
        # Create sample-exon heatmap
        pivot_z = z_scores.pivot(index='sample_id', columns='exon', values='z_score')
        
        plt.figure(figsize=(10, max(8, len(pivot_z) * 0.2)))
        
        # Custom colormap for CNV interpretation
        colors = ['darkred', 'red', 'orange', 'yellow', 'lightgreen', 'green', 
                 'lightblue', 'blue', 'darkblue']
        n_bins = 100
        cmap = plt.cm.colors.ListedColormap(
            [plt.cm.RdBu_r(i) for i in np.linspace(0, 1, n_bins)]
        )
        
        sns.heatmap(pivot_z, cmap=cmap, center=0, 
                   cbar_kws={'label': 'Z-score'}, 
                   fmt='.2f', linewidths=0.1, cbar=True)
        
        # Add threshold lines
        plt.axhline(y=0, color='black', linewidth=0.5)
        
        plt.title('Sample-Exon Z-score Heatmap')
        plt.ylabel('Sample ID')
        plt.xlabel('Exon')
        plt.xticks(rotation=45)
        plt.yticks(rotation=0, fontsize=8)
        
        plot_file = plot_dir / 'zscore_heatmap.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['zscore_heatmap'] = plot_file
    
    # Plot 3: Method comparison (if ML data available)
    if 'exon_cn' in data and any('consensus_score' in col or 'method_agreement' in col 
                                for col in data['exon_cn'].columns):
        exon_cn = data['exon_cn']
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Consensus score distribution
        if 'consensus_score' in exon_cn.columns:
            axes[0].hist(exon_cn['consensus_score'], bins=20, alpha=0.7, edgecolor='black')
            axes[0].set_title('Consensus Score Distribution')
            axes[0].set_xlabel('Consensus Score')
            axes[0].set_ylabel('Call Count')
            axes[0].axvline(x=0.8, color='red', linestyle='--', alpha=0.7, label='High consensus threshold')
            axes[0].legend()
        
        # Method agreement vs confidence
        if 'method_agreement' in exon_cn.columns and 'confidence_score' in exon_cn.columns:
            scatter = axes[1].scatter(exon_cn['method_agreement'], exon_cn['confidence_score'],
                                    alpha=0.6, c=exon_cn['copy_number'], cmap='viridis')
            axes[1].set_xlabel('Method Agreement')
            axes[1].set_ylabel('Confidence Score')
            axes[1].set_title('Method Agreement vs Confidence')
            plt.colorbar(scatter, ax=axes[1], label='Copy Number')
        
        plt.tight_layout()
        plot_file = plot_dir / 'method_comparison.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['method_comparison'] = plot_file
    
    return plots

def plot_to_base64(plot_path):
    """Convert plot to base64 for embedding in HTML."""
    try:
        with open(plot_path, 'rb') as f:
            image_data = f.read()
        return base64.b64encode(image_data).decode('utf-8')
    except Exception as e:
        print(f"Warning: Could not encode plot {plot_path}: {e}")
        return ""

def generate_html_report(data, summary, plots, output_file):
    """Generate comprehensive HTML report similar to MultiQC."""
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>SMN CNV Analysis - Consolidated Report</title>
        <style>
            body {{
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                margin: 0;
                padding: 0;
                background-color: #f5f5f5;
                color: #333;
            }}
            
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 2rem;
                text-align: center;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                padding: 2rem;
            }}
            
            .card {{
                background: white;
                border-radius: 8px;
                padding: 1.5rem;
                margin-bottom: 2rem;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            
            .card h2 {{
                margin-top: 0;
                color: #667eea;
                border-bottom: 2px solid #eee;
                padding-bottom: 0.5rem;
            }}
            
            .metrics-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 1rem;
                margin: 1rem 0;
            }}
            
            .metric-box {{
                background: #f8f9fa;
                padding: 1rem;
                border-radius: 6px;
                border-left: 4px solid #667eea;
            }}
            
            .metric-value {{
                font-size: 1.8rem;
                font-weight: bold;
                color: #667eea;
            }}
            
            .metric-label {{
                font-size: 0.9rem;
                color: #666;
                margin-top: 0.25rem;
            }}
            
            .alert {{
                padding: 1rem;
                border-radius: 6px;
                margin: 1rem 0;
            }}
            
            .alert-danger {{
                background-color: #fee;
                border-left: 4px solid #dc3545;
                color: #721c24;
            }}
            
            .alert-warning {{
                background-color: #fff3cd;
                border-left: 4px solid #ffc107;
                color: #856404;
            }}
            
            .alert-success {{
                background-color: #d1edff;
                border-left: 4px solid #28a745;
                color: #155724;
            }}
            
            .plot-container {{
                text-align: center;
                margin: 2rem 0;
            }}
            
            .plot-container img {{
                max-width: 100%;
                height: auto;
                border-radius: 6px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 1rem 0;
            }}
            
            th, td {{
                padding: 0.75rem;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }}
            
            th {{
                background-color: #f8f9fa;
                font-weight: 600;
            }}
            
            .cn-0 {{ background-color: #ffebee; color: #c62828; }}
            .cn-1 {{ background-color: #fff3e0; color: #ef6c00; }}
            .cn-2 {{ background-color: #e8f5e8; color: #2e7d32; }}
            .cn-3 {{ background-color: #e3f2fd; color: #1565c0; }}
            .cn-4 {{ background-color: #f3e5f5; color: #7b1fa2; }}
            
            .confidence-high {{ background-color: #c8e6c9; }}
            .confidence-medium {{ background-color: #ffe0b2; }}
            .confidence-low {{ background-color: #ffcdd2; }}
            
            .tabs {{
                display: flex;
                border-bottom: 2px solid #eee;
                margin-bottom: 1rem;
            }}
            
            .tab {{
                padding: 0.75rem 1.5rem;
                background: none;
                border: none;
                cursor: pointer;
                font-size: 1rem;
                color: #666;
                border-bottom: 2px solid transparent;
            }}
            
            .tab.active {{
                color: #667eea;
                border-bottom-color: #667eea;
            }}
            
            .tab-content {{
                display: none;
            }}
            
            .tab-content.active {{
                display: block;
            }}
            
            @media (max-width: 768px) {{
                .container {{
                    padding: 1rem;
                }}
                
                .metrics-grid {{
                    grid-template-columns: 1fr;
                }}
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>🧬 SMN CNV Analysis Report</h1>
            <p>Comprehensive Copy Number Variation Analysis for SMN1/SMN2 Genes</p>
            <p><small>Generated on {summary['run_info']['analysis_date']}</small></p>
        </div>
        
        <div class="container">
    """
    
    # Run Summary Card
    html_content += f"""
            <div class="card">
                <h2>📊 Run Summary</h2>
                <div class="metrics-grid">
                    <div class="metric-box">
                        <div class="metric-value">{summary['run_info']['total_samples']}</div>
                        <div class="metric-label">Total Samples</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{summary['run_info']['reference_samples']}</div>
                        <div class="metric-label">Reference Samples</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{summary['run_info']['test_samples']}</div>
                        <div class="metric-label">Test Samples</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{len(summary['run_info']['genes_analyzed'])}</div>
                        <div class="metric-label">Genes Analyzed</div>
                    </div>
                </div>
                
                <p><strong>Genes:</strong> {', '.join(summary['run_info']['genes_analyzed'])}</p>
                <p><strong>Exons:</strong> {', '.join(summary['run_info']['exons_analyzed'])}</p>
            </div>
    """
    
    # Quality Metrics Card
    if summary['quality_metrics']:
        qm = summary['quality_metrics']
        html_content += f"""
            <div class="card">
                <h2>🎯 Quality Metrics</h2>
                <div class="metrics-grid">
                    <div class="metric-box">
                        <div class="metric-value">{qm.get('mean_z_score', 0):.3f}</div>
                        <div class="metric-label">Mean Z-score</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{qm.get('std_z_score', 0):.3f}</div>
                        <div class="metric-label">Z-score Std Dev</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{qm.get('high_confidence_calls', 0)}</div>
                        <div class="metric-label">High Confidence Calls</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-value">{qm.get('low_outlier_probability_calls', 0)}</div>
                        <div class="metric-label">Low Outlier Probability</div>
                    </div>
                </div>
        """
        
        # Add quality alerts
        if abs(qm.get('mean_z_score', 0)) > 0.5:
            html_content += '<div class="alert alert-warning">⚠️ Mean Z-score deviates significantly from zero - check reference samples</div>'
        
        if qm.get('std_z_score', 0) > 2.0:
            html_content += '<div class="alert alert-warning">⚠️ High Z-score variability - check sample quality</div>'
        
        html_content += "</div>"
    
    # CNV Summary Card
    if summary['cnv_summary']:
        html_content += """
            <div class="card">
                <h2>🔬 Copy Number Variation Summary</h2>
                <div class="tabs">
                    <button class="tab active" onclick="showTab(event, 'smn1-tab')">SMN1</button>
                    <button class="tab" onclick="showTab(event, 'smn2-tab')">SMN2</button>
                </div>
        """
        
        # SMN1 Tab
        if 'SMN1' in summary['cnv_summary']:
            smn1 = summary['cnv_summary']['SMN1']
            html_content += f"""
                <div id="smn1-tab" class="tab-content active">
                    <div class="metrics-grid">
                        <div class="metric-box">
                            <div class="metric-value">{smn1.get('affected_samples', 0)}</div>
                            <div class="metric-label">Affected (CN≤0.5)</div>
                        </div>
                        <div class="metric-box">
                            <div class="metric-value">{smn1.get('carrier_samples', 0)}</div>
                            <div class="metric-label">Carriers (CN=1)</div>
                        </div>
                        <div class="metric-box">
                            <div class="metric-value">{smn1.get('normal_samples', 0)}</div>
                            <div class="metric-label">Normal (CN=2)</div>
                        </div>
                        <div class="metric-box">
                            <div class="metric-value">{smn1.get('duplication_samples', 0)

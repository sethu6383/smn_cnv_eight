#!/usr/bin/env python3

"""
generate_consolidated_report.py - Generate MultiQC-style consolidated report for SMN CNV analysis
Usage: python generate_consolidated_report.py <results_dir> <output_prefix> [--format html,txt]
"""

import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

def load_pipeline_data(results_dir):
    """Load all pipeline data files."""
    results_path = Path(results_dir)
    data = {}
    
    # Load copy number results
    cn_file = results_path / 'cnv_calls' / 'copy_numbers.txt'
    if cn_file.exists():
        data['copy_numbers'] = pd.read_csv(cn_file, sep='\t')
        print(f"Loaded copy number data: {len(data['copy_numbers'])} records")
    
    # Load gene-level results
    gene_file = results_path / 'cnv_calls' / 'copy_numbers_gene_level.txt'
    if gene_file.exists():
        data['gene_results'] = pd.read_csv(gene_file, sep='\t')
        print(f"Loaded gene-level data: {len(data['gene_results'])} records")
    
    # Load coverage data
    coverage_file = results_path / 'depth' / 'coverage_summary.txt'
    if coverage_file.exists():
        data['coverage'] = pd.read_csv(coverage_file, sep='\t')
        print(f"Loaded coverage data: {len(data['coverage'])} records")
    
    # Load Z-scores
    zscore_file = results_path / 'normalized' / 'z_scores.txt'
    if zscore_file.exists():
        data['z_scores'] = pd.read_csv(zscore_file, sep='\t')
        print(f"Loaded Z-score data: {len(data['z_scores'])} records")
    
    # Load allele counts if available
    allele_file = results_path / 'allele_counts' / 'allele_counts.txt'
    if allele_file.exists():
        data['allele_counts'] = pd.read_csv(allele_file, sep='\t')
        print(f"Loaded allele count data: {len(data['allele_counts'])} records")
    
    # Load bootstrap results if available
    bootstrap_file = results_path / 'cnv_calls' / 'copy_numbers_bootstrap.txt'
    if bootstrap_file.exists():
        data['bootstrap'] = pd.read_csv(bootstrap_file, sep='\t')
        print(f"Loaded bootstrap data: {len(data['bootstrap'])} records")
    
    # Load analysis metadata
    metadata_file = results_path / 'analysis_metadata.txt'
    if metadata_file.exists():
        with open(metadata_file, 'r') as f:
            data['metadata'] = f.read()
    
    return data

def calculate_summary_statistics(data):
    """Calculate comprehensive summary statistics."""
    stats = {
        'pipeline_info': {},
        'sample_stats': {},
        'coverage_stats': {},
        'cnv_stats': {},
        'quality_stats': {},
        'clinical_stats': {}
    }
    
    # Pipeline information
    if 'metadata' in data:
        stats['pipeline_info']['metadata'] = data['metadata']
    
    # Sample statistics
    if 'copy_numbers' in data:
        cn_df = data['copy_numbers']
        stats['sample_stats'] = {
            'total_samples': cn_df['sample_id'].nunique(),
            'total_exon_calls': len(cn_df),
            'exons_analyzed': sorted(cn_df['exon'].unique()),
            'sample_types': cn_df['sample_type'].value_counts().to_dict() if 'sample_type' in cn_df.columns else {}
        }
    
    # Coverage statistics
    if 'coverage' in data:
        cov_df = data['coverage']
        stats['coverage_stats'] = {
            'mean_coverage': cov_df['avg_coverage'].mean(),
            'median_coverage': cov_df['avg_coverage'].median(),
            'low_coverage_samples': (cov_df['avg_coverage'] < 20).sum(),
            'coverage_by_exon': cov_df.groupby('exon')['avg_coverage'].agg(['mean', 'std']).to_dict('index')
        }
    
    # CNV statistics
    if 'copy_numbers' in data:
        cn_df = data['copy_numbers']
        stats['cnv_stats'] = {
            'cn_distribution': cn_df['copy_number'].value_counts().sort_index().to_dict(),
            'confidence_distribution': cn_df['confidence'].value_counts().to_dict(),
            'category_distribution': cn_df['cn_category'].value_counts().to_dict()
        }
        
        # Per-exon CNV stats
        stats['cnv_stats']['by_exon'] = {}
        for exon in cn_df['exon'].unique():
            exon_data = cn_df[cn_df['exon'] == exon]
            stats['cnv_stats']['by_exon'][exon] = {
                'cn_distribution': exon_data['copy_number'].value_counts().sort_index().to_dict(),
                'confidence_distribution': exon_data['confidence'].value_counts().to_dict()
            }
    
    # Quality statistics
    if 'copy_numbers' in data:
        cn_df = data['copy_numbers']
        stats['quality_stats'] = {
            'high_confidence_calls': (cn_df['confidence'] == 'high').sum(),
            'low_confidence_calls': (cn_df['confidence'] == 'low').sum(),
            'consensus_score_mean': cn_df['consensus_score'].mean() if 'consensus_score' in cn_df.columns else 0,
            'method_agreement_mean': cn_df['method_agreement'].mean() if 'method_agreement' in cn_df.columns else 0
        }
        
        if 'bootstrap' in data:
            boot_df = data['bootstrap']
            stats['quality_stats']['bootstrap_confidence_mean'] = boot_df['bootstrap_confidence'].mean()
            stats['quality_stats']['low_bootstrap_confidence'] = (boot_df['bootstrap_confidence'] < 0.5).sum()
    
    # Clinical statistics (gene-level)
    if 'gene_results' in data:
        gene_df = data['gene_results']
        stats['clinical_stats'] = {}
        
        # SMN1 statistics
        smn1_data = gene_df[gene_df['gene'] == 'SMN1']
        if not smn1_data.empty:
            stats['clinical_stats']['smn1'] = {
                'total_samples': len(smn1_data),
                'clinical_significance': smn1_data['clinical_significance'].value_counts().to_dict(),
                'cn_distribution': smn1_data['estimated_copy_number'].round().value_counts().sort_index().to_dict(),
                'carrier_samples': smn1_data['clinical_significance'].str.contains('carrier', case=False, na=False).sum(),
                'affected_samples': smn1_data['clinical_significance'].str.contains('affected', case=False, na=False).sum(),
                'duplication_samples': smn1_data['clinical_significance'].str.contains('duplication', case=False, na=False).sum()
            }
        
        # SMN2 statistics
        smn2_data = gene_df[gene_df['gene'] == 'SMN2']
        if not smn2_data.empty:
            stats['clinical_stats']['smn2'] = {
                'total_samples': len(smn2_data),
                'clinical_significance': smn2_data['clinical_significance'].value_counts().to_dict(),
                'cn_distribution': smn2_data['estimated_copy_number'].round().value_counts().sort_index().to_dict()
            }
    
    return stats

def create_summary_plots(data, output_dir):
    """Create summary plots for the consolidated report."""
    plot_dir = Path(output_dir)
    plot_dir.mkdir(exist_ok=True)
    
    plots_created = []
    
    try:
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Plot 1: Copy number distribution across all samples
        if 'copy_numbers' in data:
            cn_df = data['copy_numbers']
            
            plt.figure(figsize=(12, 8))
            
            # Subplot 1: Overall CN distribution
            plt.subplot(2, 2, 1)
            cn_dist = cn_df['copy_number'].value_counts().sort_index()
            plt.bar(cn_dist.index, cn_dist.values, alpha=0.7, color='skyblue', edgecolor='black')
            plt.xlabel('Copy Number')
            plt.ylabel('Number of Calls')
            plt.title('Overall Copy Number Distribution')
            plt.xticks(range(0, int(cn_dist.index.max()) + 1))
            
            # Add count labels on bars
            for i, v in enumerate(cn_dist.values):
                plt.text(cn_dist.index[i], v + 0.5, str(v), ha='center')
            
            # Subplot 2: Confidence distribution
            plt.subplot(2, 2, 2)
            conf_dist = cn_df['confidence'].value_counts()
            colors = {'high': 'green', 'medium': 'orange', 'low': 'red'}
            bar_colors = [colors.get(conf, 'gray') for conf in conf_dist.index]
            plt.bar(conf_dist.index, conf_dist.values, color=bar_colors, alpha=0.7, edgecolor='black')
            plt.xlabel('Confidence Level')
            plt.ylabel('Number of Calls')
            plt.title('Confidence Distribution')
            plt.xticks(rotation=45)
            
            # Add percentage labels
            total_calls = len(cn_df)
            for i, (conf, count) in enumerate(conf_dist.items()):
                pct = count / total_calls * 100
                plt.text(i, count + total_calls * 0.01, f'{pct:.1f}%', ha='center')
            
            # Subplot 3: CN by exon
            plt.subplot(2, 2, 3)
            exons = sorted(cn_df['exon'].unique())
            cn_by_exon = []
            for exon in exons:
                exon_data = cn_df[cn_df['exon'] == exon]['copy_number']
                cn_by_exon.append(exon_data)
            
            plt.boxplot(cn_by_exon, labels=[e.replace('_', '\n') for e in exons])
            plt.xlabel('Exon')
            plt.ylabel('Copy Number')
            plt.title('Copy Number Distribution by Exon')
            plt.xticks(rotation=45)
            
            # Subplot 4: Coverage vs Copy Number (if coverage data available)
            plt.subplot(2, 2, 4)
            if 'coverage' in data:
                # Merge CN and coverage data
                cov_df = data['coverage']
                merged = cn_df.merge(cov_df[['sample_id', 'exon', 'avg_coverage']], 
                                   on=['sample_id', 'exon'], how='left')
                
                plt.scatter(merged['avg_coverage'], merged['copy_number'], 
                          alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
                plt.xlabel('Average Coverage')
                plt.ylabel('Copy Number')
                plt.title('Coverage vs Copy Number')
                plt.axvline(x=20, color='red', linestyle='--', alpha=0.5, label='Min coverage (20x)')
                plt.legend()
            else:
                plt.text(0.5, 0.5, 'Coverage data\nnot available', 
                        ha='center', va='center', transform=plt.gca().transAxes)
                plt.title('Coverage vs Copy Number')
            
            plt.tight_layout()
            plot_file = plot_dir / 'cnv_summary.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            plots_created.append(plot_file)
        
        # Plot 2: Clinical significance (gene-level)
        if 'gene_results' in data:
            gene_df = data['gene_results']
            
            plt.figure(figsize=(14, 6))
            
            # SMN1 clinical significance
            plt.subplot(1, 2, 1)
            smn1_data = gene_df[gene_df['gene'] == 'SMN1']
            if not smn1_data.empty:
                clin_sig = smn1_data['clinical_significance'].value_counts()
                colors = plt.cm.Set3(np.linspace(0, 1, len(clin_sig)))
                wedges, texts, autotexts = plt.pie(clin_sig.values, labels=None, autopct='%1.1f%%',
                                                  colors=colors, startangle=90)
                plt.title('SMN1 Clinical Significance')
                
                # Create legend with shortened labels
                legend_labels = [label.replace(' - ', '\n') for label in clin_sig.index]
                plt.legend(wedges, legend_labels, title="Categories", 
                          loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
            
            # SMN1 copy number distribution
            plt.subplot(1, 2, 2)
            if not smn1_data.empty:
                cn_dist = smn1_data['estimated_copy_number'].round().value_counts().sort_index()
                bars = plt.bar(cn_dist.index, cn_dist.values, alpha=0.7, 
                             color='lightcoral', edgecolor='black')
                plt.xlabel('SMN1 Copy Number')
                plt.ylabel('Number of Samples')
                plt.title('SMN1 Copy Number Distribution')
                
                # Add clinical interpretation lines
                plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='SMA affected')
                plt.axvline(x=1.5, color='orange', linestyle='--', alpha=0.7, label='SMA carrier')
                plt.legend()
                
                # Add count labels
                for bar, count in zip(bars, cn_dist.values):
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                           str(count), ha='center', va='bottom')
            
            plt.tight_layout()
            plot_file = plot_dir / 'clinical_summary.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            plots_created.append(plot_file)
        
        print(f"Created {len(plots_created)} summary plots")
        
    except Exception as e:
        print(f"Warning: Could not create all plots: {e}")
    
    return plots_created

def generate_html_report(stats, plots, output_file):
    """Generate comprehensive HTML report."""
    
    # Get plot file names for embedding
    plot_files = [p.name for p in plots] if plots else []
    
    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN CNV Analysis - Consolidated Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}
        
        .header {{
            text-align: center;
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        
        .header h1 {{
            color: #2c3e50;
            margin: 0;
            font-size: 2.5em;
        }}
        
        .header .subtitle {{
            color: #7f8c8d;
            font-size: 1.2em;
            margin: 10px 0;
        }}
        
        .section {{
            margin: 30px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
            border-left: 4px solid #3498db;
        }}
        
        .section h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .metric-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            border: 1px solid #dee2e6;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        .metric-title {{
            font-weight: bold;
            color: #495057;
            margin-bottom: 10px;
            font-size: 0.9em;
            text-transform: uppercase;
        }}
        
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 5px;
        }}
        
        .metric-subtitle {{
            color: #6c757d;
            font-size: 0.9em;
        }}
        
        .plot-container {{
            text-align: center;
            margin: 30px 0;
        }}
        
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .clinical-alert {{
            padding: 15px;
            margin: 20px 0;
            border-radius: 8px;
            border-left: 4px solid #e74c3c;
            background-color: #fdf2f2;
        }}
        
        .clinical-alert.warning {{
            border-left-color: #f39c12;
            background-color: #fef9e7;
        }}
        
        .clinical-alert.success {{
            border-left-color: #27ae60;
            background-color: #eafaf1;
        }}
        
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }}
        
        th {{
            background-color: #f8f9fa;
            font-weight: bold;
            color: #495057;
        }}
        
        .footer {{
            margin-top: 50px;
            padding-top: 20px;
            border-top: 2px solid #ecf0f1;
            text-align: center;
            color: #6c757d;
        }}
        
        .quality-indicator {{
            display: inline-block;
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin-right: 8px;
        }}
        
        .quality-high {{ background-color: #27ae60; }}
        .quality-medium {{ background-color: #f39c12; }}
        .quality-low {{ background-color: #e74c3c; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 SMN CNV Analysis</h1>
            <div class="subtitle">Consolidated Report - Enhanced Pipeline v2.0</div>
            <div class="subtitle">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
        </div>
        
        <!-- Executive Summary -->
        <div class="section">
            <h2>📊 Executive Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-title">Total Samples</div>
                    <div class="metric-value">{stats['sample_stats'].get('total_samples', 0)}</div>
                    <div class="metric-subtitle">Analyzed samples</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Exon-Level Calls</div>
                    <div class="metric-value">{stats['sample_stats'].get('total_exon_calls', 0)}</div>
                    <div class="metric-subtitle">CNV calls made</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Mean Coverage</div>
                    <div class="metric-value">{stats['coverage_stats'].get('mean_coverage', 0):.1f}x</div>
                    <div class="metric-subtitle">Across all exons</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">High Confidence</div>
                    <div class="metric-value">{stats['quality_stats'].get('high_confidence_calls', 0)}</div>
                    <div class="metric-subtitle">Reliable CNV calls</div>
                </div>
            </div>
        </div>
        
        <!-- Clinical Findings -->
        <div class="section">
            <h2>🏥 Clinical Findings (SMN1)</h2>
    """
    
    # Add SMN1 clinical findings
    if 'smn1' in stats.get('clinical_stats', {}):
        smn1 = stats['clinical_stats']['smn1']
        
        html_template += f"""
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-title">Total SMN1 Samples</div>
                    <div class="metric-value">{smn1.get('total_samples', 0)}</div>
                    <div class="metric-subtitle">Gene-level calls</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">SMA Carriers</div>
                    <div class="metric-value">{smn1.get('carrier_samples', 0)}</div>
                    <div class="metric-subtitle">Require counseling</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">SMA Affected</div>
                    <div class="metric-value">{smn1.get('affected_samples', 0)}</div>
                    <div class="metric-subtitle">Immediate attention</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Duplications</div>
                    <div class="metric-value">{smn1.get('duplication_samples', 0)}</div>
                    <div class="metric-subtitle">May be protective</div>
                </div>
            </div>
        """
        
        # Clinical alerts
        affected_count = smn1.get('affected_samples', 0)
        carrier_count = smn1.get('carrier_samples', 0)
        
        if affected_count > 0:
            html_template += f"""
            <div class="clinical-alert">
                <strong>🚨 CRITICAL:</strong> {affected_count} sample(s) with potential SMA-affecting variants detected. 
                Immediate clinical correlation and confirmation testing recommended.
            </div>
            """
        
        if carrier_count > 0:
            html_template += f"""
            <div class="clinical-alert warning">
                <strong>⚠️ IMPORTANT:</strong> {carrier_count} sample(s) identified as potential SMA carriers. 
                Genetic counseling recommended.
            </div>
            """
        
        if affected_count == 0 and carrier_count == 0:
            html_template += f"""
            <div class="clinical-alert success">
                <strong>✅ NORMAL:</strong> No SMA-affecting variants or carriers detected in this cohort.
            </div>
            """
    
    html_template += """
        </div>
        
        <!-- Quality Metrics -->
        <div class="section">
            <h2>🎯 Quality Metrics</h2>
    """
    
    # Quality metrics section
    quality_stats = stats.get('quality_stats', {})
    total_calls = stats['sample_stats'].get('total_exon_calls', 1)  # Avoid division by zero
    high_conf_pct = (quality_stats.get('high_confidence_calls', 0) / total_calls) * 100
    
    html_template += f"""
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-title">
                        <span class="quality-indicator quality-high"></span>High Confidence
                    </div>
                    <div class="metric-value">{high_conf_pct:.1f}%</div>
                    <div class="metric-subtitle">{quality_stats.get('high_confidence_calls', 0)} calls</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">
                        <span class="quality-indicator quality-low"></span>Low Confidence
                    </div>
                    <div class="metric-value">{quality_stats.get('low_confidence_calls', 0)}</div>
                    <div class="metric-subtitle">Require review</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Consensus Score</div>
                    <div class="metric-value">{quality_stats.get('consensus_score_mean', 0):.3f}</div>
                    <div class="metric-subtitle">Method agreement</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Low Coverage Samples</div>
                    <div class="metric-value">{stats['coverage_stats'].get('low_coverage_samples', 0)}</div>
                    <div class="metric-subtitle">< 20x coverage</div>
                </div>
            </div>
        </div>
    """
    
    # Add plots if available
    if plots:
        html_template += """
        <!-- Summary Plots -->
        <div class="section">
            <h2>📈 Analysis Plots</h2>
        """
        
        for plot_file in plot_files:
            html_template += f"""
            <div class="plot-container">
                <img src="{plot_file}" alt="Analysis Plot">
            </div>
            """
        
        html_template += """
        </div>
        """
    
    # Copy number distribution table
    if 'cnv_stats' in stats and 'cn_distribution' in stats['cnv_stats']:
        html_template += """
        <!-- Copy Number Distribution -->
        <div class="section">
            <h2>📊 Copy Number Distribution</h2>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Copy Number</th>
                            <th>Number of Calls</th>
                            <th>Percentage</th>
                            <th>Clinical Interpretation</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        
        cn_dist = stats['cnv_stats']['cn_distribution']
        total_calls = sum(cn_dist.values())
        
        cn_interpretation = {
            0: "Homozygous deletion",
            1: "Heterozygous deletion (SMA carrier)",
            2: "Normal copy number",
            3: "Duplication",
            4: "High amplification"
        }
        
        for cn in sorted(cn_dist.keys()):
            count = cn_dist[cn]
            pct = (count / total_calls) * 100
            interpretation = cn_interpretation.get(cn, "Unknown")
            
            html_template += f"""
                        <tr>
                            <td>{cn}</td>
                            <td>{count}</td>
                            <td>{pct:.1f}%</td>
                            <td>{interpretation}</td>
                        </tr>
            """
        
        html_template += """
                    </tbody>
                </table>
            </div>
        </div>
        """
    
    # Recommendations section
    html_template += """
        <!-- Recommendations -->
        <div class="section">
            <h2>🔍 Quality Control Recommendations</h2>
            <ul>
    """
    
    recommendations = []
    
    # Check for low confidence calls
    low_conf = quality_stats.get('low_confidence_calls', 0)
    if low_conf > 0:
        recommendations.append(f"Review {low_conf} low confidence calls manually")
    
    # Check for low coverage samples
    low_cov = stats['coverage_stats'].get('low_coverage_samples', 0)
    if low_cov > 0:
        recommendations.append(f"Investigate {low_cov} samples with coverage <20x")
    
    # Check consensus scores
    consensus_mean = quality_stats.get('consensus_score_mean', 1.0)
    if consensus_mean < 0.7:
        recommendations.append(f"Overall consensus score is {consensus_mean:.3f} - consider additional validation")
    
    # Clinical recommendations
    if 'smn1' in stats.get('clinical_stats', {}):
        smn1 = stats['clinical_stats']['smn1']
        if smn1.get('affected_samples', 0) > 0:
            recommendations.append("Confirm SMA-affecting variants with orthogonal methods")
        if smn1.get('carrier_samples', 0) > 0:
            recommendations.append("Provide genetic counseling for SMA carriers")
    
    if not recommendations:
        recommendations.append("All quality metrics are within acceptable ranges")
    
    for rec in recommendations:
        html_template += f"                <li>{rec}</li>\n"
    
    html_template += """
            </ul>
        </div>
        
        <!-- Technical Details -->
        <div class="section">
            <h2>⚙️ Technical Details</h2>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Parameter</th>
                            <th>Value</th>
                            <th>Description</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    # Add technical parameters
    technical_params = [
        ("Pipeline Version", "2.0 Enhanced (Exon 8 Focus)", "Updated pipeline focusing on SMN1/SMN2 exon 8"),
        ("Analysis Date", datetime.now().strftime('%Y-%m-%d'), "Date of analysis"),
        ("Exons Analyzed", ", ".join(stats['sample_stats'].get('exons_analyzed', [])), "Target exons"),
        ("CNV Calling Method", "Consensus (Probabilistic + Threshold)", "Multiple method consensus"),
        ("Quality Threshold", "≥Medium confidence", "Minimum reliability level"),
        ("Coverage Threshold", "≥20x", "Minimum coverage requirement")
    ]
    
    for param, value, desc in technical_params:
        html_template += f"""
                        <tr>
                            <td><strong>{param}</strong></td>
                            <td>{value}</td>
                            <td>{desc}</td>
                        </tr>
        """
    
    html_template += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>SMN CNV Detection Pipeline v2.0</strong></p>
            <p>Enhanced pipeline focusing on exon 8 for improved accuracy</p>
            <p><em>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</em></p>
            <p><small>⚠️ This report is for research purposes. Clinical decisions should be based on validated results.</small></p>
        </div>
    </div>
</body>
</html>
    """
    
    # Write HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_template)
    
    print(f"HTML report generated: {output_file}")

def generate_text_report(stats, output_file):
    """Generate comprehensive text report."""
    
    report_lines = [
        "=" * 80,
        "SMN CNV ANALYSIS - CONSOLIDATED REPORT",
        "Enhanced Pipeline v2.0 (Exon 8 Focus)",
        "=" * 80,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "EXECUTIVE SUMMARY",
        "-" * 50,
        f"Total Samples Analyzed: {stats['sample_stats'].get('total_samples', 0)}",
        f"Total Exon-Level Calls: {stats['sample_stats'].get('total_exon_calls', 0)}",
        f"Exons Analyzed: {', '.join(stats['sample_stats'].get('exons_analyzed', []))}",
        f"Mean Coverage: {stats['coverage_stats'].get('mean_coverage', 0):.1f}x",
        ""
    ]
    
    # Clinical findings section
    if 'clinical_stats' in stats and 'smn1' in stats['clinical_stats']:
        smn1 = stats['clinical_stats']['smn1']
        report_lines.extend([
            "CLINICAL FINDINGS (SMN1)",
            "-" * 50,
            f"Total SMN1 Samples: {smn1.get('total_samples', 0)}",
            f"SMA Affected (CN≤0.5): {smn1.get('affected_samples', 0)}",
            f"SMA Carriers (CN=1): {smn1.get('carrier_samples', 0)}",
            f"Gene Duplications (CN≥3): {smn1.get('duplication_samples', 0)}",
            ""
        ])
        
        # Clinical significance breakdown
        if 'clinical_significance' in smn1:
            report_lines.append("Clinical Significance Breakdown:")
            for sig, count in smn1['clinical_significance'].items():
                report_lines.append(f"  - {sig}: {count} samples")
            report_lines.append("")
    
    # Quality metrics section
    quality_stats = stats.get('quality_stats', {})
    total_calls = stats['sample_stats'].get('total_exon_calls', 1)
    
    report_lines.extend([
        "QUALITY METRICS",
        "-" * 50,
        f"High Confidence Calls: {quality_stats.get('high_confidence_calls', 0)} ({quality_stats.get('high_confidence_calls', 0)/total_calls*100:.1f}%)",
        f"Low Confidence Calls: {quality_stats.get('low_confidence_calls', 0)} ({quality_stats.get('low_confidence_calls', 0)/total_calls*100:.1f}%)",
        f"Mean Consensus Score: {quality_stats.get('consensus_score_mean', 0):.3f}",
        f"Mean Method Agreement: {quality_stats.get('method_agreement_mean', 0):.3f}",
        f"Low Coverage Samples (<20x): {stats['coverage_stats'].get('low_coverage_samples', 0)}",
        ""
    ])
    
    # Bootstrap results if available
    if 'bootstrap_confidence_mean' in quality_stats:
        report_lines.extend([
            f"Bootstrap Confidence (mean): {quality_stats['bootstrap_confidence_mean']:.3f}",
            f"Low Bootstrap Confidence: {quality_stats.get('low_bootstrap_confidence', 0)}",
            ""
        ])
    
    # Copy number distribution
    if 'cnv_stats' in stats and 'cn_distribution' in stats['cnv_stats']:
        report_lines.extend([
            "COPY NUMBER DISTRIBUTION",
            "-" * 50
        ])
        
        cn_dist = stats['cnv_stats']['cn_distribution']
        total_cn_calls = sum(cn_dist.values())
        
        cn_interpretation = {
            0: "Homozygous deletion",
            1: "Heterozygous deletion",
            2: "Normal copy number",
            3: "Duplication",
            4: "High amplification"
        }
        
        for cn in sorted(cn_dist.keys()):
            count = cn_dist[cn]
            pct = (count / total_cn_calls) * 100
            interpretation = cn_interpretation.get(cn, "Unknown")
            report_lines.append(f"CN={cn}: {count:4d} calls ({pct:5.1f}%) - {interpretation}")
        
        report_lines.append("")
    
    # Per-exon statistics
    if 'cnv_stats' in stats and 'by_exon' in stats['cnv_stats']:
        report_lines.extend([
            "PER-EXON STATISTICS",
            "-" * 50
        ])
        
        for exon, exon_stats in stats['cnv_stats']['by_exon'].items():
            report_lines.append(f"\n{exon}:")
            
            # Copy number distribution for this exon
            cn_dist = exon_stats.get('cn_distribution', {})
            if cn_dist:
                for cn in sorted(cn_dist.keys()):
                    count = cn_dist[cn]
                    report_lines.append(f"  CN={cn}: {count} calls")
            
            # Confidence distribution for this exon
            conf_dist = exon_stats.get('confidence_distribution', {})
            if conf_dist:
                report_lines.append("  Confidence:")
                for conf, count in conf_dist.items():
                    report_lines.append(f"    {conf}: {count} calls")
        
        report_lines.append("")
    
    # Coverage statistics by exon
    if 'coverage_stats' in stats and 'coverage_by_exon' in stats['coverage_stats']:
        report_lines.extend([
            "COVERAGE STATISTICS BY EXON",
            "-" * 50
        ])
        
        for exon, cov_stats in stats['coverage_stats']['coverage_by_exon'].items():
            mean_cov = cov_stats.get('mean', 0)
            std_cov = cov_stats.get('std', 0)
            report_lines.append(f"{exon}: {mean_cov:.1f}x ± {std_cov:.1f}x")
        
        report_lines.append("")
    
    # Quality control recommendations
    report_lines.extend([
        "QUALITY CONTROL RECOMMENDATIONS",
        "-" * 50
    ])
    
    recommendations = []
    
    # Check for low confidence calls
    low_conf = quality_stats.get('low_confidence_calls', 0)
    if low_conf > 0:
        recommendations.append(f"• Review {low_conf} low confidence calls manually")
    
    # Check for low coverage samples
    low_cov = stats['coverage_stats'].get('low_coverage_samples', 0)
    if low_cov > 0:
        recommendations.append(f"• Investigate {low_cov} samples with coverage <20x")
    
    # Check consensus scores
    consensus_mean = quality_stats.get('consensus_score_mean', 1.0)
    if consensus_mean < 0.7:
        recommendations.append(f"• Overall consensus score is {consensus_mean:.3f} - consider additional validation")
    
    # Clinical recommendations
    if 'smn1' in stats.get('clinical_stats', {}):
        smn1 = stats['clinical_stats']['smn1']
        if smn1.get('affected_samples', 0) > 0:
            recommendations.append("• Confirm SMA-affecting variants with orthogonal methods")
        if smn1.get('carrier_samples', 0) > 0:
            recommendations.append("• Provide genetic counseling for SMA carriers")
    
    if not recommendations:
        recommendations.append("• All quality metrics are within acceptable ranges")
    
    report_lines.extend(recommendations)
    report_lines.extend(["", ""])
    
    # Technical details
    report_lines.extend([
        "TECHNICAL DETAILS",
        "-" * 50,
        f"Pipeline Version: 2.0 Enhanced (Exon 8 Focus)",
        f"Analysis Date: {datetime.now().strftime('%Y-%m-%d')}",
        f"CNV Calling Method: Consensus (Probabilistic + Threshold)",
        f"Quality Threshold: ≥Medium confidence",
        f"Coverage Threshold: ≥20x",
        "",
        "IMPORTANT NOTES",
        "-" * 50,
        "• This pipeline focuses exclusively on exon 8 for improved accuracy",
        "• All CNV calls should be validated through orthogonal methods", 
        "• Clinical interpretation requires correlation with patient phenotype",
        "• Consider family studies for carrier samples when appropriate",
        "",
        "=" * 80,
        f"Report generated by SMN CNV Detection Pipeline v2.0",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "⚠️  This report is for research purposes. Clinical decisions should be based on validated results.",
        "=" * 80
    ])
    
    # Write text file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Text report generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate consolidated MultiQC-style report for SMN CNV analysis')
    parser.add_argument('results_dir', help='Results directory containing pipeline outputs')
    parser.add_argument('output_prefix', help='Output file prefix (e.g., /path/to/SMN_CNV_Analysis)')
    parser.add_argument('--format', default='html,txt', help='Output formats: html, txt, or both (default: html,txt)')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading pipeline data...")
    data = load_pipeline_data(args.results_dir)
    
    if not data:
        print("Error: No pipeline data found in results directory")
        sys.exit(1)
    
    print("Calculating summary statistics...")
    stats = calculate_summary_statistics(data)
    
    # Create plots
    print("Creating summary plots...")
    plot_output_dir = output_dir / 'plots'
    plots = create_summary_plots(data, plot_output_dir)
    
    # Generate reports based on format
    formats = args.format.lower().split(',')
    
    if 'html' in formats:
        html_file = f"{args.output_prefix}_consolidated_report.html"
        print(f"Generating HTML report: {html_file}")
        generate_html_report(stats, plots, html_file)
    
    if 'txt' in formats:
        txt_file = f"{args.output_prefix}_consolidated_report.txt"
        print(f"Generating text report: {txt_file}")
        generate_text_report(stats, txt_file)
    
    print("\nConsolidated report generation completed successfully!")
    
    # Print summary to console
    print(f"\nSUMMARY:")
    print(f"  Total samples: {stats['sample_stats'].get('total_samples', 0)}")
    print(f"  Total CNV calls: {stats['sample_stats'].get('total_exon_calls', 0)}")
    print(f"  Mean coverage: {stats['coverage_stats'].get('mean_coverage', 0):.1f}x")
    
    if 'smn1' in stats.get('clinical_stats', {}):
        smn1 = stats['clinical_stats']['smn1']
        print(f"  SMN1 affected: {smn1.get('affected_samples', 0)}")
        print(f"  SMN1 carriers: {smn1.get('carrier_samples', 0)}")
    
    high_conf = stats['quality_stats'].get('high_confidence_calls', 0)
    total_calls = stats['sample_stats'].get('total_exon_calls', 1)
    print(f"  High confidence: {high_conf}/{total_calls} ({high_conf/total_calls*100:.1f}%)")

if __name__ == "__main__":
    main()

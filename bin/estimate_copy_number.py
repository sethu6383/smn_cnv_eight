#!/usr/bin/env python3

"""
estimate_copy_number.py - Enhanced CNV estimation with ML features and bootstrap confidence
Usage: python estimate_copy_number.py <z_scores_file> <output_file> [--thresholds custom_thresholds.txt] [--bootstrap-samples 1000]
"""

import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix
import joblib
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Enhanced copy number thresholds with confidence intervals
DEFAULT_THRESHOLDS = {
    'homozygous_deletion': -2.5,    # CN=0
    'heterozygous_deletion': -1.5,  # CN=1
    'normal_lower': -1.5,           # CN=2 lower bound
    'normal_upper': 1.5,            # CN=2 upper bound
    'duplication': 2.5,             # CN=3
    'high_amplification': 3.5,      # CN=4+
}

# ML-enhanced thresholds (can be learned from data)
ML_THRESHOLDS = {
    'confidence_threshold': 0.7,    # Minimum confidence for high-confidence calls
    'outlier_threshold': 0.1,       # Maximum outlier probability for reliable calls
    'consensus_threshold': 0.8,     # Minimum consensus score for final calls
}

def read_custom_thresholds(threshold_file):
    """Read custom thresholds from file."""
    thresholds = DEFAULT_THRESHOLDS.copy()
    
    try:
        with open(threshold_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    threshold_name, value = parts[0], float(parts[1])
                    if threshold_name in thresholds:
                        thresholds[threshold_name] = value
        
        print(f"Loaded custom thresholds from: {threshold_file}")
    except Exception as e:
        print(f"Warning: Could not read custom thresholds ({e}). Using defaults.")
    
    return thresholds

def bootstrap_z_score_thresholds(z_scores_df, reference_samples, n_bootstrap=1000):
    """Use bootstrap sampling to estimate confidence intervals for thresholds."""
    ref_data = z_scores_df[z_scores_df['sample_id'].isin(reference_samples)]
    
    if len(ref_data) < 10:
        print("Warning: Insufficient reference samples for bootstrap analysis")
        return DEFAULT_THRESHOLDS
    
    bootstrap_stats = []
    
    for _ in range(n_bootstrap):
        # Bootstrap sample
        boot_sample = ref_data.sample(n=len(ref_data), replace=True)
        
        # Calculate statistics for this bootstrap sample
        for exon in boot_sample['exon'].unique():
            exon_data = boot_sample[boot_sample['exon'] == exon]['z_score']
            if len(exon_data) > 5:
                stats_dict = {
                    'exon': exon,
                    'mean': exon_data.mean(),
                    'std': exon_data.std(),
                    'q05': exon_data.quantile(0.05),
                    'q95': exon_data.quantile(0.95),
                    'q025': exon_data.quantile(0.025),
                    'q975': exon_data.quantile(0.975)
                }
                bootstrap_stats.append(stats_dict)
    
    if bootstrap_stats:
        bootstrap_df = pd.DataFrame(bootstrap_stats)
        
        # Calculate confidence intervals for thresholds
        enhanced_thresholds = DEFAULT_THRESHOLDS.copy()
        
        # Use bootstrap quantiles to adjust thresholds
        for exon in bootstrap_df['exon'].unique():
            exon_bootstrap = bootstrap_df[bootstrap_df['exon'] == exon]
            
            # Calculate exon-specific thresholds based on bootstrap CI
            q025_mean = exon_bootstrap['q025'].mean()
            q975_mean = exon_bootstrap['q975'].mean()
            
            # Adjust thresholds based on reference distribution
            enhanced_thresholds[f'{exon}_deletion_threshold'] = q025_mean
            enhanced_thresholds[f'{exon}_duplication_threshold'] = q975_mean
        
        return enhanced_thresholds
    
    return DEFAULT_THRESHOLDS

def probabilistic_cnv_calling(z_score, robust_z_score, confidence_score, outlier_prob, 
                             thresholds, ml_predictions=None):
    """Enhanced probabilistic CNV calling using multiple evidence sources."""
    
    # Base probability from Z-score
    z_prob = calculate_z_score_probability(z_score, thresholds)
    
    # Robust Z-score probability
    robust_prob = calculate_z_score_probability(robust_z_score, thresholds)
    
    # Confidence weighting
    confidence_weight = min(1.0, max(0.1, confidence_score))
    
    # Outlier penalty
    outlier_penalty = max(0.0, 1.0 - outlier_prob * 5)  # Penalize high outlier probability
    
    # Combine probabilities
    base_prob = (z_prob + robust_prob) / 2
    weighted_prob = base_prob * confidence_weight * outlier_penalty
    
    # ML enhancement if available
    if ml_predictions:
        ml_weight = 0.3  # Weight for ML predictions
        if 'random_forest_max_prob' in ml_predictions:
            ml_confidence = ml_predictions['random_forest_max_prob']
            weighted_prob = weighted_prob * (1 - ml_weight) + ml_confidence * ml_weight
    
    # Convert probability to copy number
    cn_estimate, cn_category, final_confidence = probability_to_copy_number(
        weighted_prob, z_score, thresholds
    )
    
    return {
        'copy_number': cn_estimate,
        'cn_category': cn_category,
        'confidence': final_confidence,
        'base_probability': base_prob,
        'weighted_probability': weighted_prob,
        'confidence_weight': confidence_weight,
        'outlier_penalty': outlier_penalty
    }

def calculate_z_score_probability(z_score, thresholds):
    """Calculate copy number probabilities from Z-score."""
    if pd.isna(z_score):
        return {'CN0': 0.25, 'CN1': 0.25, 'CN2': 0.25, 'CN3+': 0.25}
    
    # Define probability distributions for each copy number state
    # Using normal distributions centered at expected Z-scores
    cn_distributions = {
        'CN0': {'mean': -3.0, 'std': 0.5},    # Homozygous deletion
        'CN1': {'mean': -2.0, 'std': 0.3},    # Heterozygous deletion  
        'CN2': {'mean': 0.0, 'std': 0.5},     # Normal
        'CN3+': {'mean': 2.5, 'std': 0.7}     # Duplication/amplification
    }
    
    probabilities = {}
    for cn, dist in cn_distributions.items():
        prob = stats.norm.pdf(z_score, dist['mean'], dist['std'])
        probabilities[cn] = prob
    
    # Normalize probabilities
    total_prob = sum(probabilities.values())
    if total_prob > 0:
        probabilities = {k: v/total_prob for k, v in probabilities.items()}
    
    return probabilities

def probability_to_copy_number(prob_dict, z_score, thresholds):
    """Convert probability distribution to discrete copy number call."""
    if isinstance(prob_dict, dict):
        # Find most likely copy number
        max_cn = max(prob_dict.keys(), key=lambda k: prob_dict[k])
        max_prob = prob_dict[max_cn]
        
        # Map to integer copy number
        cn_mapping = {'CN0': 0, 'CN1': 1, 'CN2': 2, 'CN3+': 3}
        cn_estimate = cn_mapping.get(max_cn, 2)
        
        # Determine category
        category_mapping = {
            0: 'homozygous_deletion',
            1: 'heterozygous_deletion', 
            2: 'normal',
            3: 'duplication'
        }
        cn_category = category_mapping.get(cn_estimate, 'unknown')
        
        # Confidence based on probability and consistency
        if max_prob > 0.8:
            confidence = 'high'
        elif max_prob > 0.6:
            confidence = 'medium'
        else:
            confidence = 'low'
        
    else:
        # Fallback to threshold-based calling
        if pd.isna(z_score):
            return np.nan, 'unknown', 'low'
        
        if z_score <= thresholds['homozygous_deletion']:
            cn_estimate, cn_category = 0, 'homozygous_deletion'
        elif z_score <= thresholds['heterozygous_deletion']:
            cn_estimate, cn_category = 1, 'heterozygous_deletion'
        elif thresholds['normal_lower'] < z_score <= thresholds['normal_upper']:
            cn_estimate, cn_category = 2, 'normal'
        elif z_score <= thresholds['duplication']:
            cn_estimate, cn_category = 3, 'duplication'
        else:
            cn_estimate, cn_category = 4, 'high_amplification'
        
        # Distance-based confidence
        boundary_distances = [
            abs(z_score - thresholds['homozygous_deletion']),
            abs(z_score - thresholds['heterozygous_deletion']),
            abs(z_score - thresholds['normal_lower']),
            abs(z_score - thresholds['normal_upper']),
            abs(z_score - thresholds['duplication'])
        ]
        
        min_distance = min(boundary_distances)
        if min_distance > 1.0:
            confidence = 'high'
        elif min_distance > 0.5:
            confidence = 'medium'
        else:
            confidence = 'low'
    
    return cn_estimate, cn_category, confidence

def consensus_cnv_calling(row, thresholds, ml_predictions_available=False):
    """Consensus CNV calling using multiple methods."""
    methods = []
    
    # Method 1: Standard Z-score
    std_result = probabilistic_cnv_calling(
        row['z_score'], row.get('z_score_robust', row['z_score']),
        row.get('confidence_score', 0.5), row.get('outlier_probability', 0.5),
        thresholds
    )
    methods.append(('standard', std_result))
    
    # Method 2: Robust Z-score (if available)
    if 'z_score_robust' in row:
        robust_result = probabilistic_cnv_calling(
            row['z_score_robust'], row['z_score_robust'],
            row.get('confidence_score', 0.5), row.get('outlier_probability', 0.5),
            thresholds
        )
        methods.append(('robust', robust_result))
    
    # Method 3: ML predictions (if available)
    if ml_predictions_available and 'random_forest_prediction' in row:
        ml_result = {
            'copy_number': row['random_forest_prediction'],
            'cn_category': get_category_from_cn(row['random_forest_prediction']),
            'confidence': 'high' if row.get('random_forest_max_prob', 0) > 0.8 else 'medium',
            'base_probability': row.get('random_forest_max_prob', 0.5),
            'weighted_probability': row.get('random_forest_max_prob', 0.5),
            'confidence_weight': 1.0,
            'outlier_penalty': 1.0
        }
        methods.append(('ml', ml_result))
    
    # Consensus logic
    cn_votes = [method[1]['copy_number'] for method in methods if not pd.isna(method[1]['copy_number'])]
    
    if not cn_votes:
        return {
            'copy_number': np.nan,
            'cn_category': 'unknown',
            'confidence': 'low',
            'consensus_score': 0.0,
            'method_agreement': 0.0,
            'contributing_methods': []
        }
    
    # Calculate consensus
    from collections import Counter
    cn_counts = Counter(cn_votes)
    consensus_cn = cn_counts.most_common(1)[0][0]
    consensus_score = cn_counts[consensus_cn] / len(cn_votes)
    
    # Average confidence from contributing methods
    confidences = [method[1]['confidence'] for method in methods if method[1]['copy_number'] == consensus_cn]
    conf_scores = {'high': 1.0, 'medium': 0.6, 'low': 0.3}
    avg_confidence_score = np.mean([conf_scores.get(conf, 0.3) for conf in confidences])
    
    if avg_confidence_score > 0.8:
        final_confidence = 'high'
    elif avg_confidence_score > 0.5:
        final_confidence = 'medium'
    else:
        final_confidence = 'low'
    
    return {
        'copy_number': consensus_cn,
        'cn_category': get_category_from_cn(consensus_cn),
        'confidence': final_confidence,
        'consensus_score': consensus_score,
        'method_agreement': consensus_score,
        'contributing_methods': [method[0] for method in methods if method[1]['copy_number'] == consensus_cn]
    }

def get_category_from_cn(cn):
    """Map copy number to category."""
    if pd.isna(cn):
        return 'unknown'
    elif cn == 0:
        return 'homozygous_deletion'
    elif cn == 1:
        return 'heterozygous_deletion'
    elif cn == 2:
        return 'normal'
    elif cn == 3:
        return 'duplication'
    else:
        return 'high_amplification'

def bootstrap_confidence_intervals(z_scores_df, n_bootstrap=1000):
    """Calculate bootstrap confidence intervals for copy number estimates."""
    bootstrap_results = []
    
    samples = z_scores_df['sample_id'].unique()
    
    for sample_id in samples:
        sample_data = z_scores_df[z_scores_df['sample_id'] == sample_id]
        
        sample_bootstrap = []
        for _ in range(n_bootstrap):
            # Bootstrap resample (with replacement) of the sample's data points
            boot_sample = sample_data.sample(n=len(sample_data), replace=True)
            
            # Calculate mean Z-score for this bootstrap sample
            mean_z = boot_sample['z_score'].mean()
            sample_bootstrap.append(mean_z)
        
        # Calculate confidence intervals
        ci_lower = np.percentile(sample_bootstrap, 2.5)
        ci_upper = np.percentile(sample_bootstrap, 97.5)
        ci_width = ci_upper - ci_lower
        
        bootstrap_results.append({
            'sample_id': sample_id,
            'bootstrap_mean': np.mean(sample_bootstrap),
            'bootstrap_std': np.std(sample_bootstrap),
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'ci_width': ci_width,
            'bootstrap_confidence': 1.0 / (1.0 + ci_width)  # Narrow CI = high confidence
        })
    
    return pd.DataFrame(bootstrap_results)

def estimate_gene_copy_numbers_enhanced(cn_results_df):
    """Enhanced gene-level copy number estimation with uncertainty quantification."""
    gene_results = []
    
    for sample_id in cn_results_df['sample_id'].unique():
        sample_data = cn_results_df[cn_results_df['sample_id'] == sample_id]
        
        # Separate SMN1 and SMN2 based on exon names  
        smn1_data = sample_data[sample_data['exon'].str.contains('SMN1')]
        smn2_data = sample_data[sample_data['exon'].str.contains('SMN2')]
        
        for gene_name, gene_data in [('SMN1', smn1_data), ('SMN2', smn2_data)]:
            if not gene_data.empty:
                # Multiple estimation methods
                cn_values = gene_data['copy_number'].dropna()
                consensus_scores = gene_data['consensus_score'].fillna(0.5)
                method_agreements = gene_data['method_agreement'].fillna(0.5)
                
                if len(cn_values) > 0:
                    # Weighted average based on consensus scores
                    weights = consensus_scores.values
                    weighted_cn = np.average(cn_values, weights=weights)
                    
                    # Uncertainty metrics
                    cn_std = cn_values.std() if len(cn_values) > 1 else 0
                    mean_consensus = consensus_scores.mean()
                    mean_agreement = method_agreements.mean()
                    
                    # Final confidence assessment
                    confidence_factors = [
                        mean_consensus,           # How well methods agree
                        mean_agreement,           # Method agreement score
                        1.0 / (1.0 + cn_std),   # Consistency across exons
                        len(cn_values) / 4.0     # Number of supporting exons
                    ]
                    
                    final_confidence_score = np.mean(confidence_factors)
                    
                    if final_confidence_score > 0.8:
                        final_confidence = 'high'
                    elif final_confidence_score > 0.6:
                        final_confidence = 'medium'
                    else:
                        final_confidence = 'low'
                    
                    # Clinical significance
                    clinical_significance = assess_clinical_significance(
                        gene_name, weighted_cn, final_confidence_score
                    )
                    
                    gene_results.append({
                        'sample_id': sample_id,
                        'gene': gene_name,
                        'estimated_copy_number': round(weighted_cn, 2),
                        'copy_number_std': round(cn_std, 3),
                        'consensus_score': round(mean_consensus, 3),
                        'method_agreement': round(mean_agreement, 3),
                        'confidence': final_confidence,
                        'confidence_score': round(final_confidence_score, 3),
                        'clinical_significance': clinical_significance,
                        'n_exons': len(gene_data),
                        'supporting_exons': len(cn_values)
                    })
    
    return pd.DataFrame(gene_results)

def assess_clinical_significance(gene, copy_number, confidence_score):
    """Assess clinical significance of copy number variant."""
    if gene == 'SMN1':
        if copy_number <= 0.5 and confidence_score > 0.7:
            return 'Pathogenic - SMA affected'
        elif copy_number <= 1.5 and confidence_score > 0.6:
            return 'Likely pathogenic - SMA carrier'
        elif 1.5 < copy_number < 2.5:
            return 'Benign - Normal copy number'
        elif copy_number >= 2.5:
            return 'Uncertain - Gene duplication'
        else:
            return 'Uncertain - Low confidence'
    
    elif gene == 'SMN2':
        if copy_number == 0 and confidence_score > 0.7:
            return 'Modifier - SMN2 deletion may worsen SMA'
        elif copy_number >= 3 and confidence_score > 0.6:
            return 'Modifier - SMN2 duplication may improve SMA'
        else:
            return 'Benign - Normal SMN2 variation'
    
    return 'Unknown'

def create_enhanced_visualization(cn_results_df, gene_results_df, output_dir, thresholds, bootstrap_results=None):
    """Create comprehensive visualization suite."""
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Enhanced copy number distribution with confidence
    plt.figure(figsize=(15, 10))
    
    exons = sorted(cn_results_df['exon'].unique())
    
    for i, exon in enumerate(exons):
        plt.subplot(2, 3, i+1)
        exon_data = cn_results_df[cn_results_df['exon'] == exon]
        
        # Color by confidence
        high_conf = exon_data[exon_data['confidence'] == 'high']
        med_conf = exon_data[exon_data['confidence'] == 'medium']
        low_conf = exon_data[exon_data['confidence'] == 'low']
        
        plt.hist(high_conf['copy_number'], alpha=0.7, label='High confidence', color='green', bins=range(0, 6))
        plt.hist(med_conf['copy_number'], alpha=0.7, label='Medium confidence', color='orange', bins=range(0, 6))
        plt.hist(low_conf['copy_number'], alpha=0.7, label='Low confidence', color='red', bins=range(0, 6))
        
        plt.title(f'{exon} Copy Number Distribution')
        plt.xlabel('Copy Number')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.xticks(range(0, 5))
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'enhanced_cn_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Consensus scores vs copy number
    plt.figure(figsize=(12, 8))
    
    for i, exon in enumerate(exons):
        plt.subplot(1, len(exons), i+1)
        exon_data = cn_results_df[cn_results_df['exon'] == exon]
        
        scatter = plt.scatter(exon_data['copy_number'], exon_data['consensus_score'],
                            c=exon_data['method_agreement'], cmap='viridis',
                            alpha=0.7, s=50)
        plt.xlabel('Copy Number')
        plt.ylabel('Consensus Score')
        plt.title(f'{exon}')
        plt.colorbar(scatter, label='Method Agreement')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'consensus_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Gene-level clinical significance
    if not gene_results_df.empty:
        plt.figure(figsize=(12, 8))
        
        # Create significance category mapping
        sig_categories = gene_results_df['clinical_significance'].value_counts()
        
        plt.subplot(1, 2, 1)
        plt.pie(sig_categories.values, labels=sig_categories.index, autopct='%1.1f%%')
        plt.title('Clinical Significance Distribution')
        
        plt.subplot(1, 2, 2)
        
        # Confidence vs copy number for SMN1
        smn1_data = gene_results_df[gene_results_df['gene'] == 'SMN1']
        if not smn1_data.empty:
            scatter = plt.scatter(smn1_data['estimated_copy_number'], 
                                smn1_data['confidence_score'],
                                c=pd.Categorical(smn1_data['clinical_significance']).codes,
                                cmap='Set1', alpha=0.7, s=100)
            plt.xlabel('SMN1 Copy Number')
            plt.ylabel('Confidence Score')
            plt.title('SMN1: Copy Number vs Confidence')
            
            # Add clinical thresholds
            plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='SMA affected')
            plt.axvline(x=1.5, color='orange', linestyle='--', alpha=0.5, label='SMA carrier')
            plt.legend()
        
        plt.tight_layout()
        plt.savefig(plot_dir / 'clinical_significance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Plot 4: Bootstrap confidence intervals (if available)
    if bootstrap_results is not None and not bootstrap_results.empty:
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.hist(bootstrap_results['bootstrap_confidence'], bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Bootstrap Confidence')
        plt.ylabel('Sample Count')
        plt.title('Bootstrap Confidence Distribution')
        
        plt.subplot(1, 2, 2)
        plt.scatter(bootstrap_results['ci_width'], bootstrap_results['bootstrap_confidence'])
        plt.xlabel('Confidence Interval Width')
        plt.ylabel('Bootstrap Confidence')
        plt.title('CI Width vs Bootstrap Confidence')
        
        plt.tight_layout()
        plt.savefig(plot_dir / 'bootstrap_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Enhanced visualizations saved to: {plot_dir}")

def main():
    parser = argparse.ArgumentParser(description='Enhanced CNV estimation with ML and bootstrap confidence')
    parser.add_argument('z_scores_file', help='Z-scores file from normalize_coverage.py')
    parser.add_argument('output_file', help='Output file for copy number estimates')
    parser.add_argument('--thresholds', help='Custom thresholds file', default=None)
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    parser.add_argument('--bootstrap-samples', type=int, default=1000, help='Number of bootstrap samples')
    parser.add_argument('--enable-ml', action='store_true', help='Enable ML-enhanced predictions')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load thresholds
    if args.thresholds:
        thresholds = read_custom_thresholds(args.thresholds)
    else:
        thresholds = DEFAULT_THRESHOLDS.copy()
        print("Using default thresholds")
    
    print(f"Enhanced copy number thresholds:")
    for name, value in thresholds.items():
        print(f"  {name}: {value}")
    
    # Read Z-scores
    try:
        z_scores_df = pd.read_csv(args.z_scores_file, sep='\t')
        print(f"Loaded Z-scores for {len(z_scores_df)} sample-exon combinations")
    except Exception as e:
        print(f"Error reading Z-scores file: {e}")
        sys.exit(1)
    
    # Identify reference samples for bootstrap
    reference_samples = z_scores_df[z_scores_df['sample_type'] == 'reference']['sample_id'].unique()
    
    # Bootstrap threshold optimization (if enough reference samples)
    if len(reference_samples) >= 10:
        print("Optimizing thresholds using bootstrap analysis...")
        bootstrap_thresholds = bootstrap_z_score_thresholds(z_scores_df, reference_samples, args.bootstrap_samples)
        thresholds.update(bootstrap_thresholds)
    
    # Check for ML predictions
    ml_columns = [col for col in z_scores_df.columns if 'random_forest' in col or 'gmm' in col]
    ml_predictions_available = len(ml_columns) > 0 and args.enable_ml
    
    print("Performing enhanced CNV estimation...")
    cn_results = []
    
    for _, row in z_scores_df.iterrows():
        # Enhanced consensus calling
        cn_result = consensus_cnv_calling(row, thresholds, ml_predictions_available)
        
        result = {
            'sample_id': row['sample_id'],
            'exon': row['exon'],
            'z_score': row['z_score'],
            'z_score_robust': row.get('z_score_robust', row['z_score']),
            'copy_number': cn_result['copy_number'],
            'cn_category': cn_result['cn_category'],
            'confidence': cn_result['confidence'],
            'consensus_score': cn_result.get('consensus_score', 0.5),
            'method_agreement': cn_result.get('method_agreement', 0.5),
            'contributing_methods': ','.join(cn_result.get('contributing_methods', [])),
            'raw_coverage': row.get('raw_coverage', np.nan),
            'confidence_score': row.get('confidence_score', 0.5),
            'outlier_probability': row.get('outlier_probability', 0.5),
            'sample_type': row.get('sample_type', 'unknown')
        }
        
        # Add ML predictions if available
        if ml_predictions_available:
            for col in ml_columns:
                if col in row:
                    result[col] = row[col]
        
        cn_results.append(result)
    
    cn_results_df = pd.DataFrame(cn_results)
    
    # Bootstrap confidence intervals
    bootstrap_results = None
    if args.bootstrap_samples > 0:
        print(f"Calculating bootstrap confidence intervals ({args.bootstrap_samples} samples)...")
        bootstrap_results = bootstrap_confidence_intervals(z_scores_df, args.bootstrap_samples)
        
        # Merge bootstrap results
        cn_results_df = cn_results_df.merge(
            bootstrap_results[['sample_id', 'bootstrap_confidence', 'ci_width']], 
            on='sample_id', how='left'
        )
    
    # Enhanced gene-level estimates
    print("Calculating enhanced gene-level copy numbers...")
    gene_results_df = estimate_gene_copy_numbers_enhanced(cn_results_df)
    
    # Save results
    cn_results_df.to_csv(args.output_file, index=False, sep='\t')
    
    gene_output_file = args.output_file.replace('.txt', '_gene_level.txt')
    gene_results_df.to_csv(gene_output_file, index=False, sep='\t')
    
    # Save enhanced thresholds used
    threshold_output_file = args.output_file.replace('.txt', '_thresholds.txt')
    with open(threshold_output_file, 'w') as f:
        f.write("# Enhanced copy number thresholds used\n")
        f.write("threshold_name\tvalue\tdescription\n")
        for name, value in thresholds.items():
            f.write(f"{name}\t{value}\tEnhanced threshold\n")
    
    # Save bootstrap results if available
    if bootstrap_results is not None:
        bootstrap_output_file = args.output_file.replace('.txt', '_bootstrap.txt')
        bootstrap_results.to_csv(bootstrap_output_file, index=False, sep='\t')
    
    # Create enhanced visualizations
    if not args.no_plots:
        try:
            create_enhanced_visualization(cn_results_df, gene_results_df, output_dir, 
                                        thresholds, bootstrap_results)
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Print enhanced summary
    print(f"\nEnhanced CNV estimation completed!")
    print(f"Exon-level results saved to: {args.output_file}")
    print(f"Gene-level results saved to: {gene_output_file}")
    print(f"Thresholds saved to: {threshold_output_file}")
    
    if bootstrap_results is not None:
        print(f"Bootstrap results saved to: {bootstrap_output_file}")
    
    # Enhanced summary statistics
    print(f"\nEnhanced summary statistics:")
    print(f"Total samples: {len(cn_results_df['sample_id'].unique())}")
    
    # Confidence distribution
    conf_dist = cn_results_df['confidence'].value_counts()
    print(f"\nConfidence distribution:")
    for conf, count in conf_dist.items():
        print(f"  {conf}: {count} calls ({count/len(cn_results_df)*100:.1f}%)")
    
    # Method agreement statistics
    if 'method_agreement' in cn_results_df.columns:
        mean_agreement = cn_results_df['method_agreement'].mean()
        print(f"\nMean method agreement: {mean_agreement:.3f}")
    
    # Clinical significance summary (SMN1 focus)
    if not gene_results_df.empty:
        smn1_results = gene_results_df[gene_results_df['gene'] == 'SMN1']
        if not smn1_results.empty:
            print(f"\nSMN1 Clinical Significance:")
            clin_sig = smn1_results['clinical_significance'].value_counts()
            for sig, count in clin_sig.items():
                print(f"  {sig}: {count} samples")
    
    # Quality control recommendations
    print(f"\nQuality Control Recommendations:")
    low_conf_samples = cn_results_df[cn_results_df['confidence'] == 'low']['sample_id'].nunique()
    if low_conf_samples > 0:
        print(f"  - {low_conf_samples} samples have low confidence calls - consider manual review")
    
    if 'consensus_score' in cn_results_df.columns:
        low_consensus = (cn_results_df['consensus_score'] < 0.7).sum()
        if low_consensus > 0:
            print(f"  - {low_consensus} calls have low consensus scores - verify results")
    
    if bootstrap_results is not None:
        low_bootstrap_conf = (bootstrap_results['bootstrap_confidence'] < 0.5).sum()
        if low_bootstrap_conf > 0:
            print(f"  - {low_bootstrap_conf} samples have low bootstrap confidence")
    
    print(f"\nRecommended thresholds for reliable CNV calling:")
    print(f"  - Minimum consensus score: 0.7")
    print(f"  - Minimum confidence: medium or high")
    print(f"  - Maximum outlier probability: 0.1")
    if bootstrap_results is not None:
        print(f"  - Minimum bootstrap confidence: 0.5")

if __name__ == "__main__":
    main()

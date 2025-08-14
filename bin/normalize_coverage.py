#!/usr/bin/env python3

"""
normalize_coverage.py - Enhanced coverage normalization with ML features for exon 8 CNV detection
Usage: python normalize_coverage.py <coverage_file> <sample_info_file> <output_file> [--enable-ml]
"""

import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report
import joblib
import warnings
warnings.filterwarnings('ignore')

def read_sample_info(sample_info_file):
    """Read sample information from file."""
    try:
        sample_df = pd.read_csv(sample_info_file, sep='\t')
        samples = {}
        for _, row in sample_df.iterrows():
            samples[row['sample_id']] = {
                'bam_path': row.get('bam_path', ''),
                'sample_type': row.get('sample_type', 'unknown')
            }
        return samples
    except Exception as e:
        print(f"Warning: Could not read sample info file ({e}). Using defaults.")
        return {}

def robust_outlier_detection(data, method='iqr', contamination=0.1):
    """Detect outliers using multiple robust methods."""
    data_array = np.array(data).reshape(-1, 1)
    
    if method == 'iqr':
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        outliers = (data < lower_bound) | (data > upper_bound)
    
    elif method == 'isolation_forest':
        iso_forest = IsolationForest(contamination=contamination, random_state=42)
        outliers = iso_forest.fit_predict(data_array) == -1
    
    elif method == 'zscore':
        z_scores = np.abs(stats.zscore(data))
        outliers = z_scores > 3
    
    elif method == 'modified_zscore':
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        modified_z_scores = 0.6745 * (data - median) / mad
        outliers = np.abs(modified_z_scores) > 3.5
    
    else:
        outliers = np.zeros(len(data), dtype=bool)
    
    return outliers

def calculate_robust_reference_stats(coverage_df, reference_samples):
    """Calculate robust statistics for each exon using reference samples with advanced outlier detection."""
    ref_stats = {}
    
    # Filter for reference samples only
    ref_df = coverage_df[coverage_df['sample_id'].isin(reference_samples)]
    
    if ref_df.empty:
        print("Warning: No reference samples found in coverage data!")
        return ref_stats
    
    # Calculate statistics for each exon
    exons = ref_df['exon'].unique()
    
    for exon in exons:
        exon_data = ref_df[ref_df['exon'] == exon]['avg_coverage'].values
        
        if len(exon_data) > 0:
            # Multi-method outlier detection
            outlier_methods = ['iqr', 'isolation_forest', 'modified_zscore']
            outlier_consensus = np.zeros(len(exon_data), dtype=bool)
            
            for method in outlier_methods:
                try:
                    outliers = robust_outlier_detection(exon_data, method=method)
                    outlier_consensus |= outliers
                except:
                    continue
            
            # Use consensus outlier detection (majority vote)
            outlier_votes = sum([robust_outlier_detection(exon_data, method=method) 
                               for method in outlier_methods])
            consensus_outliers = outlier_votes >= 2  # At least 2 methods agree
            
            # Filter outliers
            filtered_data = exon_data[~consensus_outliers]
            
            if len(filtered_data) >= 3:  # Need at least 3 samples for robust stats
                # Robust statistics
                robust_mean = np.mean(filtered_data)
                robust_std = np.std(filtered_data, ddof=1)
                robust_median = np.median(filtered_data)
                
                # Additional robust estimators
                trimmed_mean = stats.trim_mean(filtered_data, 0.1)  # 10% trimmed mean
                mad = stats.median_abs_deviation(filtered_data)  # Median Absolute Deviation
                
                ref_stats[exon] = {
                    'mean': robust_mean,
                    'std': robust_std,
                    'median': robust_median,
                    'trimmed_mean': trimmed_mean,
                    'mad': mad,
                    'n_samples': len(filtered_data),
                    'n_outliers': len(exon_data) - len(filtered_data),
                    'min': np.min(filtered_data),
                    'max': np.max(filtered_data),
                    'q25': np.percentile(filtered_data, 25),
                    'q75': np.percentile(filtered_data, 75),
                    'cv': robust_std / robust_mean if robust_mean > 0 else np.inf  # Coefficient of variation
                }
            else:
                # Fallback to basic stats if too few samples
                ref_stats[exon] = {
                    'mean': np.mean(exon_data),
                    'std': np.std(exon_data, ddof=1) if len(exon_data) > 1 else 1.0,
                    'median': np.median(exon_data),
                    'trimmed_mean': np.mean(exon_data),
                    'mad': stats.median_abs_deviation(exon_data) if len(exon_data) > 1 else 1.0,
                    'n_samples': len(exon_data),
                    'n_outliers': 0,
                    'min': np.min(exon_data),
                    'max': np.max(exon_data),
                    'q25': np.percentile(exon_data, 25),
                    'q75': np.percentile(exon_data, 75),
                    'cv': np.std(exon_data) / np.mean(exon_data) if np.mean(exon_data) > 0 else np.inf
                }
        else:
            print(f"Warning: No coverage data for exon {exon} in reference samples")
    
    return ref_stats

def calculate_enhanced_z_scores(coverage_df, ref_stats, allele_df=None):
    """Calculate enhanced Z-scores with multiple normalization methods."""
    z_score_data = []
    
    # Initialize scalers
    robust_scaler = RobustScaler()
    standard_scaler = StandardScaler()
    
    for _, row in coverage_df.iterrows():
        sample_id = row['sample_id']
        exon = row['exon']
        coverage = row['avg_coverage']
        
        if exon in ref_stats:
            ref_stats_exon = ref_stats[exon]
            ref_mean = ref_stats_exon['mean']
            ref_std = ref_stats_exon['std']
            ref_median = ref_stats_exon['median']
            ref_mad = ref_stats_exon['mad']
            trimmed_mean = ref_stats_exon['trimmed_mean']
            
            # Multiple Z-score calculations
            # Standard Z-score
            z_score_standard = (coverage - ref_mean) / ref_std if ref_std > 0 else 0.0
            
            # Robust Z-score using MAD
            z_score_robust = (coverage - ref_median) / (1.4826 * ref_mad) if ref_mad > 0 else 0.0
            
            # Trimmed mean Z-score
            z_score_trimmed = (coverage - trimmed_mean) / ref_std if ref_std > 0 else 0.0
            
            # Log-transformed Z-score (for non-negative data)
            if coverage > 0 and ref_mean > 0:
                log_coverage = np.log1p(coverage)
                log_ref_mean = np.log1p(ref_mean)
                log_ref_std = np.log1p(ref_std)
                z_score_log = (log_coverage - log_ref_mean) / log_ref_std if log_ref_std > 0 else 0.0
            else:
                z_score_log = 0.0
            
            # Confidence score based on reference sample size and CV
            confidence_score = min(1.0, ref_stats_exon['n_samples'] / 10.0) * (1.0 / (1.0 + ref_stats_exon['cv']))
            
            # Outlier probability
            outlier_prob = calculate_outlier_probability(coverage, ref_stats_exon)
            
            result = {
                'sample_id': sample_id,
                'exon': exon,
                'raw_coverage': coverage,
                'ref_mean': ref_mean,
                'ref_std': ref_std,
                'ref_median': ref_median,
                'ref_mad': ref_mad,
                'z_score': z_score_standard,  # Primary Z-score
                'z_score_robust': z_score_robust,
                'z_score_trimmed': z_score_trimmed,
                'z_score_log': z_score_log,
                'ref_n_samples': ref_stats_exon['n_samples'],
                'confidence_score': confidence_score,
                'outlier_probability': outlier_prob,
                'ref_cv': ref_stats_exon['cv']
            }
            
            # Add allele fraction features if available
            if allele_df is not None and not allele_df.empty:
                sample_alleles = allele_df[allele_df['sample_id'] == sample_id]
                if not sample_alleles.empty:
                    # Calculate allele fraction features
                    allele_features = calculate_allele_features(sample_alleles, exon)
                    result.update(allele_features)
            
            z_score_data.append(result)
        else:
            print(f"Warning: No reference statistics for exon {exon}")
    
    return pd.DataFrame(z_score_data)

def calculate_outlier_probability(coverage, ref_stats):
    """Calculate probability that a coverage value is an outlier."""
    mean = ref_stats['mean']
    std = ref_stats['std']
    
    if std <= 0:
        return 0.0
    
    # Use normal distribution to calculate probability
    z_score = abs(coverage - mean) / std
    
    # Two-tailed probability
    prob = 2 * (1 - stats.norm.cdf(z_score))
    
    return min(1.0, prob)

def calculate_allele_features(allele_data, exon):
    """Calculate allele fraction features for ML models."""
    features = {
        'allele_smn1_fraction': 0.0,
        'allele_smn2_fraction': 0.0,
        'allele_total_depth': 0,
        'allele_imbalance': 0.0,
        'allele_quality_score': 0.0
    }
    
    if allele_data.empty:
        return features
    
    # Calculate features from allele data
    smn1_data = allele_data[allele_data['gene'] == 'SMN1']
    smn2_data = allele_data[allele_data['gene'] == 'SMN2']
    
    if not smn1_data.empty:
        features['allele_smn1_fraction'] = smn1_data['alt_freq'].mean()
        features['allele_total_depth'] += smn1_data['total_depth'].sum()
    
    if not smn2_data.empty:
        features['allele_smn2_fraction'] = smn2_data['alt_freq'].mean()
        features['allele_total_depth'] += smn2_data['total_depth'].sum()
    
    # Calculate allele imbalance
    if features['allele_smn1_fraction'] + features['allele_smn2_fraction'] > 0:
        features['allele_imbalance'] = abs(features['allele_smn1_fraction'] - features['allele_smn2_fraction'])
    
    # Quality score based on depth
    features['allele_quality_score'] = min(1.0, features['allele_total_depth'] / 100.0)
    
    return features

def train_ml_models(z_scores_df, reference_samples):
    """Train ML models for enhanced CNV prediction."""
    models = {}
    
    # Prepare training data
    ref_data = z_scores_df[z_scores_df['sample_id'].isin(reference_samples)]
    
    if len(ref_data) < 5:
        print("Warning: Insufficient reference samples for ML training")
        return models
    
    # Feature columns for ML
    feature_cols = ['z_score', 'z_score_robust', 'z_score_trimmed', 'z_score_log', 
                   'confidence_score', 'outlier_probability', 'ref_cv']
    
    # Add allele features if available
    allele_cols = ['allele_smn1_fraction', 'allele_smn2_fraction', 'allele_imbalance', 'allele_quality_score']
    available_allele_cols = [col for col in allele_cols if col in z_scores_df.columns]
    feature_cols.extend(available_allele_cols)
    
    # Filter available features
    available_features = [col for col in feature_cols if col in z_scores_df.columns]
    
    if len(available_features) < 3:
        print("Warning: Insufficient features for ML training")
        return models
    
    X_ref = ref_data[available_features].fillna(0)
    
    # 1. Gaussian Mixture Model for copy number clustering
    try:
        gmm = GaussianMixture(n_components=3, random_state=42)  # CN=0,1,2+ clusters
        gmm.fit(X_ref)
        models['gmm'] = {
            'model': gmm,
            'features': available_features,
            'type': 'unsupervised'
        }
        print("Trained Gaussian Mixture Model")
    except Exception as e:
        print(f"Warning: GMM training failed: {e}")
    
    # 2. Isolation Forest for anomaly detection
    try:
        iso_forest = IsolationForest(contamination=0.1, random_state=42)
        iso_forest.fit(X_ref)
        models['isolation_forest'] = {
            'model': iso_forest,
            'features': available_features,
            'type': 'anomaly_detection'
        }
        print("Trained Isolation Forest")
    except Exception as e:
        print(f"Warning: Isolation Forest training failed: {e}")
    
    # 3. If we have labeled data (from known CNVs), train supervised models
    # For now, we'll create synthetic labels based on Z-scores for demonstration
    try:
        # Create synthetic labels based on Z-score thresholds
        y_synthetic = np.zeros(len(X_ref))
        y_synthetic[ref_data['z_score'] <= -2.5] = 0  # CN=0
        y_synthetic[(ref_data['z_score'] > -2.5) & (ref_data['z_score'] <= -1.5)] = 1  # CN=1
        y_synthetic[(ref_data['z_score'] > -1.5) & (ref_data['z_score'] <= 1.5)] = 2  # CN=2
        y_synthetic[ref_data['z_score'] > 1.5] = 3  # CN=3+
        
        # Only train if we have multiple classes
        if len(np.unique(y_synthetic)) > 1:
            rf_model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
            rf_model.fit(X_ref, y_synthetic)
            
            # Calculate cross-validation score
            cv_scores = cross_val_score(rf_model, X_ref, y_synthetic, cv=min(5, len(X_ref)))
            
            models['random_forest'] = {
                'model': rf_model,
                'features': available_features,
                'type': 'supervised',
                'cv_score': cv_scores.mean(),
                'feature_importance': dict(zip(available_features, rf_model.feature_importances_))
            }
            print(f"Trained Random Forest (CV score: {cv_scores.mean():.3f})")
    except Exception as e:
        print(f"Warning: Random Forest training failed: {e}")
    
    return models

def apply_ml_predictions(z_scores_df, models):
    """Apply trained ML models to make enhanced CNV predictions."""
    if not models:
        return z_scores_df
    
    for model_name, model_info in models.items():
        try:
            model = model_info['model']
            features = model_info['features']
            model_type = model_info['type']
            
            # Prepare data
            X = z_scores_df[features].fillna(0)
            
            if model_type == 'unsupervised':
                if model_name == 'gmm':
                    predictions = model.predict(X)
                    probabilities = model.predict_proba(X)
                    z_scores_df[f'{model_name}_cluster'] = predictions
                    z_scores_df[f'{model_name}_max_prob'] = np.max(probabilities, axis=1)
                
            elif model_type == 'anomaly_detection':
                anomaly_scores = model.decision_function(X)
                anomaly_labels = model.predict(X)
                z_scores_df[f'{model_name}_anomaly_score'] = anomaly_scores
                z_scores_df[f'{model_name}_is_anomaly'] = (anomaly_labels == -1)
                
            elif model_type == 'supervised':
                predictions = model.predict(X)
                probabilities = model.predict_proba(X)
                z_scores_df[f'{model_name}_prediction'] = predictions
                z_scores_df[f'{model_name}_max_prob'] = np.max(probabilities, axis=1)
                
                # Add class probabilities
                for i, class_label in enumerate(model.classes_):
                    z_scores_df[f'{model_name}_prob_cn{class_label}'] = probabilities[:, i]
        
        except Exception as e:
            print(f"Warning: Could not apply {model_name} predictions: {e}")
    
    return z_scores_df

def create_enhanced_plots(z_scores_df, ref_stats, output_dir, models=None):
    """Create enhanced visualization plots."""
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Enhanced Z-score distributions
    plt.figure(figsize=(15, 10))
    exons = z_scores_df['exon'].unique()
    
    score_types = ['z_score', 'z_score_robust', 'z_score_trimmed']
    available_scores = [score for score in score_types if score in z_scores_df.columns]
    
    for i, exon in enumerate(sorted(exons)):
        exon_data = z_scores_df[z_scores_df['exon'] == exon]
        
        for j, score_type in enumerate(available_scores):
            plt.subplot(len(exons), len(available_scores), i * len(available_scores) + j + 1)
            
            scores = exon_data[score_type].dropna()
            if len(scores) > 0:
                plt.hist(scores, bins=20, alpha=0.7, edgecolor='black')
                plt.axvline(x=0, color='red', linestyle='--', alpha=0.7)
                plt.axvline(x=-1.5, color='orange', linestyle='--', alpha=0.5)
                plt.axvline(x=1.5, color='orange', linestyle='--', alpha=0.5)
                plt.title(f'{exon} - {score_type}')
                plt.xlabel('Score')
                plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'enhanced_score_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Confidence vs Z-score scatter
    plt.figure(figsize=(12, 8))
    
    for i, exon in enumerate(sorted(exons)):
        plt.subplot(1, len(exons), i + 1)
        exon_data = z_scores_df[z_scores_df['exon'] == exon]
        
        scatter = plt.scatter(exon_data['z_score'], exon_data['confidence_score'], 
                            c=exon_data['outlier_probability'], cmap='viridis',
                            alpha=0.7, s=50)
        plt.xlabel('Z-score')
        plt.ylabel('Confidence Score')
        plt.title(f'{exon}')
        plt.colorbar(scatter, label='Outlier Probability')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'confidence_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: ML model results (if available)
    if models and 'random_forest' in models:
        rf_info = models['random_forest']
        feature_importance = rf_info.get('feature_importance', {})
        
        if feature_importance:
            plt.figure(figsize=(10, 6))
            features = list(feature_importance.keys())
            importances = list(feature_importance.values())
            
            plt.barh(features, importances)
            plt.xlabel('Feature Importance')
            plt.title('Random Forest Feature Importance')
            plt.tight_layout()
            plt.savefig(plot_dir / 'feature_importance.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    print(f"Enhanced plots saved to: {plot_dir}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python normalize_coverage.py <coverage_file> <sample_info_file> <output_file> [--enable-ml] [--allele-file allele_file]")
        sys.exit(1)
    
    coverage_file = sys.argv[1]
    sample_info_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Parse optional arguments
    enable_ml = '--enable-ml' in sys.argv
    allele_file = None
    if '--allele-file' in sys.argv:
        allele_idx = sys.argv.index('--allele-file')
        if allele_idx + 1 < len(sys.argv):
            allele_file = sys.argv[allele_idx + 1]
    
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Reading input files...")
    
    # Read coverage data
    try:
        coverage_df = pd.read_csv(coverage_file, sep='\t')
        print(f"Loaded coverage data for {len(coverage_df)} sample-exon combinations")
    except Exception as e:
        print(f"Error reading coverage file: {e}")
        sys.exit(1)
    
    # Read allele data if available
    allele_df = None
    if allele_file and Path(allele_file).exists():
        try:
            allele_df = pd.read_csv(allele_file, sep='\t')
            print(f"Loaded allele data for {len(allele_df)} records")
        except Exception as e:
            print(f"Warning: Could not read allele file: {e}")
    
    # Read sample information
    samples_info = read_sample_info(sample_info_file)
    
    # Auto-detect if no sample info
    if not samples_info:
        print("Auto-detecting sample types from sample names...")
        for sample_id in coverage_df['sample_id'].unique():
            if any(keyword in sample_id.lower() for keyword in ['ref', 'control', 'normal']):
                sample_type = 'reference'
            else:
                sample_type = 'test'
            samples_info[sample_id] = {'sample_type': sample_type, 'bam_path': ''}
    
    # Identify reference samples
    reference_samples = [sid for sid, info in samples_info.items() 
                        if info['sample_type'] == 'reference']
    
    print(f"Found {len(reference_samples)} reference samples")
    
    if len(reference_samples) < 3:
        print("Warning: Very few reference samples. Results may be unreliable.")
    
    # Calculate enhanced reference statistics
    print("Calculating robust reference statistics...")
    ref_stats = calculate_robust_reference_stats(coverage_df, reference_samples)
    
    # Calculate enhanced Z-scores
    print("Calculating enhanced Z-scores...")
    z_scores_df = calculate_enhanced_z_scores(coverage_df, ref_stats, allele_df)
    
    # Add sample type information
    z_scores_df['sample_type'] = z_scores_df['sample_id'].map(
        lambda x: samples_info.get(x, {}).get('sample_type', 'unknown')
    )
    
    # Apply ML models if enabled
    models = {}
    if enable_ml:
        print("Training ML models...")
        models = train_ml_models(z_scores_df, reference_samples)
        
        if models:
            print("Applying ML predictions...")
            z_scores_df = apply_ml_predictions(z_scores_df, models)
            
            # Save models
            models_dir = output_dir / 'models'
            models_dir.mkdir(exist_ok=True)
            
            for model_name, model_info in models.items():
                try:
                    model_file = models_dir / f'{model_name}_model.joblib'
                    joblib.dump(model_info, model_file)
                    print(f"Saved {model_name} model to {model_file}")
                except Exception as e:
                    print(f"Warning: Could not save {model_name} model: {e}")
    
    # Save results
    z_scores_df.to_csv(output_file, index=False, sep='\t')
    
    # Save enhanced reference statistics
    ref_stats_file = output_file.replace('.txt', '_ref_stats.txt')
    ref_stats_df = pd.DataFrame.from_dict(ref_stats, orient='index').reset_index()
    ref_stats_df.rename(columns={'index': 'exon'}, inplace=True)
    ref_stats_df.to_csv(ref_stats_file, index=False, sep='\t')
    
    # Create enhanced plots
    try:
        create_enhanced_plots(z_scores_df, ref_stats, output_dir, models)
    except Exception as e:
        print(f"Warning: Could not create plots: {e}")
    
    # Print enhanced summary
    print(f"\nEnhanced normalization completed!")
    print(f"Z-scores saved to: {output_file}")
    print(f"Reference statistics saved to: {ref_stats_file}")
    
    if enable_ml and models:
        print(f"ML models trained: {list(models.keys())}")
        print(f"ML models saved to: {output_dir}/models/")
    
    print(f"\nEnhanced summary statistics:")
    print(f"Total samples: {len(z_scores_df['sample_id'].unique())}")
    print(f"Reference samples: {len(reference_samples)}")
    
    # Enhanced statistics
    for exon in z_scores_df['exon'].unique():
        exon_data = z_scores_df[z_scores_df['exon'] == exon]
        print(f"\n{exon} statistics:")
        print(f"  Samples: {len(exon_data)}")
        print(f"  Z-score range: {exon_data['z_score'].min():.2f} to {exon_data['z_score'].max():.2f}")
        print(f"  Mean confidence: {exon_data['confidence_score'].mean():.3f}")
        print(f"  High outlier probability (>0.1): {(exon_data['outlier_probability'] > 0.1).sum()}")

if __name__ == "__main__":
    main()

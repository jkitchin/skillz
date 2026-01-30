"""
Comprehensive Outlier Detection Example

This example demonstrates:
- Multiple outlier detection methods (statistical, PyOD, sklearn)
- Ensemble detection (combining multiple methods)
- Visualization of outlier scores
- Investigation and decision-making for outliers
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# PyOD detectors
from pyod.models.iforest import IForest
from pyod.models.lof import LOF
from pyod.models.copod import COPOD
from pyod.models.ecod import ECOD
from pyod.models.knn import KNN
from pyod.models.pca import PCA as PCA_OD

# sklearn detectors
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.preprocessing import StandardScaler


def generate_data_with_outliers(n_inliers=500, n_outliers=25, seed=42):
    """Generate sample data with known outliers"""
    np.random.seed(seed)

    # Generate inliers (normal data)
    X_inliers = np.random.randn(n_inliers, 2) * [2, 1] + [5, 3]

    # Generate outliers
    X_outliers = np.random.uniform(low=-5, high=15, size=(n_outliers, 2))

    # Combine
    X = np.vstack([X_inliers, X_outliers])
    true_labels = np.hstack([np.zeros(n_inliers), np.ones(n_outliers)])

    return X, true_labels


def detect_outliers_statistical(X):
    """Univariate statistical outlier detection methods"""
    print("=" * 80)
    print("STATISTICAL OUTLIER DETECTION")
    print("=" * 80)

    results = {}

    for i, col_name in enumerate(["Feature 1", "Feature 2"]):
        data = X[:, i]

        # Z-score method
        z_scores = np.abs(stats.zscore(data))
        outliers_z = z_scores > 3

        # Modified Z-score (MAD)
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        modified_z = 0.6745 * (data - median) / mad
        outliers_mad = np.abs(modified_z) > 3.5

        # IQR method
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        outliers_iqr = (data < lower) | (data > upper)

        print(f"\n{col_name}:")
        print(f"  Z-score method:       {outliers_z.sum()} outliers")
        print(f"  Modified Z-score (MAD): {outliers_mad.sum()} outliers")
        print(f"  IQR method:           {outliers_iqr.sum()} outliers")

        results[col_name] = {"z_score": outliers_z, "mad": outliers_mad, "iqr": outliers_iqr}

    # Combine: point is outlier if flagged in either feature
    combined_z = results["Feature 1"]["z_score"] | results["Feature 2"]["z_score"]
    combined_mad = results["Feature 1"]["mad"] | results["Feature 2"]["mad"]
    combined_iqr = results["Feature 1"]["iqr"] | results["Feature 2"]["iqr"]

    print(f"\nCombined (either feature):")
    print(f"  Z-score:       {combined_z.sum()} outliers")
    print(f"  Modified Z-score: {combined_mad.sum()} outliers")
    print(f"  IQR:           {combined_iqr.sum()} outliers")

    return {"z_score": combined_z, "mad": combined_mad, "iqr": combined_iqr}


def detect_outliers_pyod(X, contamination=0.1):
    """PyOD outlier detection methods"""
    print("\n" + "=" * 80)
    print(f"PyOD OUTLIER DETECTION (contamination={contamination})")
    print("=" * 80)

    # Scale data for distance-based methods
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Define detectors
    detectors = {
        "Isolation Forest": IForest(contamination=contamination, random_state=42),
        "LOF": LOF(contamination=contamination, n_neighbors=20),
        "KNN": KNN(contamination=contamination, n_neighbors=5),
        "COPOD": COPOD(contamination=contamination),
        "ECOD": ECOD(contamination=contamination),
        "PCA": PCA_OD(contamination=contamination),
    }

    results = {}
    scores_dict = {}

    print()
    for name, detector in detectors.items():
        detector.fit(X_scaled)
        labels = detector.labels_  # 0=inlier, 1=outlier
        scores = detector.decision_scores_

        results[name] = labels
        scores_dict[name] = scores

        print(f"{name:20s}: {labels.sum():3d} outliers detected")

    return results, scores_dict


def detect_outliers_sklearn(X, contamination=0.1):
    """scikit-learn outlier detection methods"""
    print("\n" + "=" * 80)
    print(f"SKLEARN OUTLIER DETECTION (contamination={contamination})")
    print("=" * 80)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    results = {}

    # Isolation Forest
    iso = IsolationForest(contamination=contamination, random_state=42)
    labels_iso = iso.fit_predict(X_scaled)  # -1=outlier, 1=inlier
    results["IsolationForest"] = (labels_iso == -1).astype(int)

    # LOF (novelty detection mode for consistency)
    lof = LocalOutlierFactor(contamination=contamination, novelty=True)
    lof.fit(X_scaled)
    labels_lof = lof.predict(X_scaled)
    results["LOF"] = (labels_lof == -1).astype(int)

    print()
    for name, labels in results.items():
        print(f"{name:20s}: {labels.sum():3d} outliers detected")

    return results


def ensemble_detection(all_results):
    """Combine multiple detection methods"""
    print("\n" + "=" * 80)
    print("ENSEMBLE OUTLIER DETECTION")
    print("=" * 80)

    # Convert all to numpy arrays
    labels_array = np.array([labels for labels in all_results.values()]).T

    # Count how many methods flagged each point
    vote_counts = labels_array.sum(axis=1)

    # Different consensus thresholds
    consensus_any = vote_counts >= 1  # At least 1 method
    consensus_2 = vote_counts >= 2  # At least 2 methods
    consensus_3 = vote_counts >= 3  # At least 3 methods
    consensus_majority = vote_counts >= len(all_results) / 2

    print(f"\nTotal methods: {len(all_results)}")
    print(f"\nConsensus results:")
    print(f"  Any method (≥1):        {consensus_any.sum()} outliers")
    print(f"  At least 2 methods:     {consensus_2.sum()} outliers")
    print(f"  At least 3 methods:     {consensus_3.sum()} outliers")
    print(f"  Majority (≥{len(all_results) // 2 + 1}):       {consensus_majority.sum()} outliers")

    return {
        "vote_counts": vote_counts,
        "consensus_any": consensus_any,
        "consensus_2": consensus_2,
        "consensus_3": consensus_3,
        "consensus_majority": consensus_majority,
    }


def evaluate_detection(detected_labels, true_labels):
    """Evaluate outlier detection performance"""
    print("\n" + "=" * 80)
    print("DETECTION PERFORMANCE EVALUATION")
    print("=" * 80)

    from sklearn.metrics import (
        confusion_matrix,
        precision_score,
        recall_score,
        f1_score,
        accuracy_score,
    )

    # Confusion matrix
    cm = confusion_matrix(true_labels, detected_labels)
    tn, fp, fn, tp = cm.ravel()

    # Metrics
    precision = precision_score(true_labels, detected_labels)
    recall = recall_score(true_labels, detected_labels)
    f1 = f1_score(true_labels, detected_labels)
    accuracy = accuracy_score(true_labels, detected_labels)

    print("\nConfusion Matrix:")
    print(f"  True Negatives:  {tn}")
    print(f"  False Positives: {fp}")
    print(f"  False Negatives: {fn}")
    print(f"  True Positives:  {tp}")

    print("\nMetrics:")
    print(f"  Precision: {precision:.3f} (of flagged points, how many are actual outliers)")
    print(f"  Recall:    {recall:.3f} (of actual outliers, how many were detected)")
    print(f"  F1-Score:  {f1:.3f} (harmonic mean of precision and recall)")
    print(f"  Accuracy:  {accuracy:.3f}")

    return {"precision": precision, "recall": recall, "f1": f1, "accuracy": accuracy}


def plot_outlier_scores(X, scores_dict, true_labels):
    """Plot outlier score distributions"""
    n_detectors = len(scores_dict)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()

    for i, (name, scores) in enumerate(scores_dict.items()):
        ax = axes[i]

        # Histogram for inliers and outliers
        inlier_scores = scores[true_labels == 0]
        outlier_scores = scores[true_labels == 1]

        ax.hist(inlier_scores, bins=30, alpha=0.6, label="Inliers", color="blue")
        ax.hist(outlier_scores, bins=30, alpha=0.6, label="Outliers", color="red")

        ax.set_xlabel("Outlier Score")
        ax.set_ylabel("Frequency")
        ax.set_title(name)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_detection_results(X, true_labels, detected_labels, method_name):
    """Visualize detection results in 2D"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # True labels
    ax = axes[0]
    ax.scatter(
        X[true_labels == 0, 0],
        X[true_labels == 0, 1],
        c="blue",
        s=50,
        alpha=0.6,
        label="True Inliers",
    )
    ax.scatter(
        X[true_labels == 1, 0],
        X[true_labels == 1, 1],
        c="red",
        s=100,
        alpha=0.8,
        marker="x",
        label="True Outliers",
    )
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.set_title("True Labels")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Detected labels
    ax = axes[1]
    ax.scatter(
        X[detected_labels == 0, 0],
        X[detected_labels == 0, 1],
        c="blue",
        s=50,
        alpha=0.6,
        label="Detected Inliers",
    )
    ax.scatter(
        X[detected_labels == 1, 0],
        X[detected_labels == 1, 1],
        c="red",
        s=100,
        alpha=0.8,
        marker="x",
        label="Detected Outliers",
    )
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.set_title(f"Detected Labels ({method_name})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_vote_distribution(vote_counts, true_labels, n_methods):
    """Plot distribution of votes"""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Separate by true labels
    votes_inliers = vote_counts[true_labels == 0]
    votes_outliers = vote_counts[true_labels == 1]

    x = np.arange(n_methods + 1)
    width = 0.35

    # Count votes
    inlier_counts = [np.sum(votes_inliers == i) for i in range(n_methods + 1)]
    outlier_counts = [np.sum(votes_outliers == i) for i in range(n_methods + 1)]

    ax.bar(x - width / 2, inlier_counts, width, label="True Inliers", alpha=0.7)
    ax.bar(x + width / 2, outlier_counts, width, label="True Outliers", alpha=0.7)

    ax.set_xlabel("Number of Methods Flagging Point as Outlier")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Votes Across Methods")
    ax.set_xticks(x)
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    return fig


def main():
    """Main execution function"""
    print("Generating data with known outliers...")
    X, true_labels = generate_data_with_outliers(n_inliers=500, n_outliers=25)

    print(f"Dataset: {X.shape[0]} points")
    print(f"  Inliers:  {(true_labels == 0).sum()}")
    print(f"  Outliers: {(true_labels == 1).sum()}")
    print(f"  True contamination: {(true_labels == 1).sum() / len(true_labels):.2%}\n")

    # Statistical methods
    stat_results = detect_outliers_statistical(X)

    # PyOD methods
    pyod_results, scores_dict = detect_outliers_pyod(X, contamination=0.05)

    # sklearn methods
    sklearn_results = detect_outliers_sklearn(X, contamination=0.05)

    # Combine all results
    all_results = {**stat_results, **pyod_results, **sklearn_results}

    # Ensemble detection
    ensemble_results = ensemble_detection(all_results)

    # Evaluate performance of different consensus strategies
    print("\n" + "=" * 80)
    print("PERFORMANCE COMPARISON")
    print("=" * 80)

    strategies = {
        "Consensus ≥2": ensemble_results["consensus_2"],
        "Consensus ≥3": ensemble_results["consensus_3"],
        "Majority Vote": ensemble_results["consensus_majority"],
        "IForest (PyOD)": pyod_results["Isolation Forest"],
        "LOF (PyOD)": pyod_results["LOF"],
    }

    performance = []
    for name, labels in strategies.items():
        metrics = evaluate_detection(labels, true_labels)
        performance.append(
            {
                "Strategy": name,
                "Precision": metrics["precision"],
                "Recall": metrics["recall"],
                "F1-Score": metrics["f1"],
            }
        )

    perf_df = pd.DataFrame(performance)
    print("\n", perf_df.to_string(index=False))

    # Visualizations
    print("\n" + "=" * 80)
    print("CREATING VISUALIZATIONS")
    print("=" * 80)

    # Outlier score distributions
    fig1 = plot_outlier_scores(X, scores_dict, true_labels)

    # Detection results (using consensus ≥2)
    fig2 = plot_detection_results(
        X, true_labels, ensemble_results["consensus_2"], "Consensus ≥2 Methods"
    )

    # Vote distribution
    fig3 = plot_vote_distribution(ensemble_results["vote_counts"], true_labels, len(all_results))

    # Show all plots
    print("\nCOMPLETE - Displaying all plots")
    plt.show()


if __name__ == "__main__":
    main()

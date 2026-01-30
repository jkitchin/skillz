"""
Robust Regression Methods for Handling Outliers

Demonstrates:
- OLS vs Robust regression (RLM, Huber, RANSAC, Theil-Sen)
- Impact of outliers on different methods
- When to use each robust method
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.robust.robust_linear_model import RLM
from sklearn.linear_model import HuberRegressor, RANSACRegressor, TheilSenRegressor
from sklearn.preprocessing import StandardScaler


def generate_data_with_outliers(n=100, n_outliers=10, seed=42):
    """Generate regression data with outliers"""
    np.random.seed(seed)

    # Inliers: y = 3 + 2*x + noise
    x_inliers = np.random.uniform(0, 10, n - n_outliers)
    y_inliers = 3 + 2 * x_inliers + np.random.normal(0, 1, n - n_outliers)

    # Outliers: random y values
    x_outliers = np.random.uniform(0, 10, n_outliers)
    y_outliers = np.random.uniform(-5, 25, n_outliers)

    # Combine
    x = np.concatenate([x_inliers, x_outliers])
    y = np.concatenate([y_inliers, y_outliers])
    is_outlier = np.concatenate([np.zeros(n - n_outliers), np.ones(n_outliers)])

    return x, y, is_outlier


def compare_regression_methods(x, y):
    """Compare OLS with robust regression methods"""
    print("=" * 80)
    print("COMPARING REGRESSION METHODS")
    print("=" * 80)

    X = sm.add_constant(x)
    results = {}

    # 1. OLS (Ordinary Least Squares)
    ols = sm.OLS(y, X).fit()
    results["OLS"] = ols
    print(f"\nOLS:")
    print(f"  Intercept: {ols.params[0]:.3f}")
    print(f"  Slope:     {ols.params[1]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean(ols.resid**2)):.3f}")

    # 2. RLM with Huber's T
    rlm_huber = RLM(y, X, M=sm.robust.norms.HuberT()).fit()
    results["RLM (Huber)"] = rlm_huber
    print(f"\nRLM (Huber's T):")
    print(f"  Intercept: {rlm_huber.params[0]:.3f}")
    print(f"  Slope:     {rlm_huber.params[1]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean(rlm_huber.resid**2)):.3f}")
    print(f"  Points with weight < 0.5: {(rlm_huber.weights < 0.5).sum()}")

    # 3. RLM with Tukey Biweight
    rlm_tukey = RLM(y, X, M=sm.robust.norms.TukeyBiweight()).fit()
    results["RLM (Tukey)"] = rlm_tukey
    print(f"\nRLM (Tukey Biweight):")
    print(f"  Intercept: {rlm_tukey.params[0]:.3f}")
    print(f"  Slope:     {rlm_tukey.params[1]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean(rlm_tukey.resid**2)):.3f}")
    print(f"  Points with weight < 0.5: {(rlm_tukey.weights < 0.5).sum()}")

    # 4. Huber Regressor (sklearn)
    huber = HuberRegressor(epsilon=1.35, max_iter=200)
    huber.fit(x.reshape(-1, 1), y)
    y_pred_huber = huber.predict(x.reshape(-1, 1))
    results["Huber (sklearn)"] = huber
    print(f"\nHuber Regressor (sklearn):")
    print(f"  Intercept: {huber.intercept_:.3f}")
    print(f"  Slope:     {huber.coef_[0]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean((y - y_pred_huber) ** 2)):.3f}")

    # 5. RANSAC (very robust, finds inlier subset)
    ransac = RANSACRegressor(min_samples=0.5, residual_threshold=2.0, random_state=42)
    ransac.fit(x.reshape(-1, 1), y)
    y_pred_ransac = ransac.predict(x.reshape(-1, 1))
    results["RANSAC"] = ransac
    inlier_mask = ransac.inlier_mask_
    print(f"\nRANSAC:")
    print(f"  Intercept: {ransac.estimator_.intercept_:.3f}")
    print(f"  Slope:     {ransac.estimator_.coef_[0]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean((y - y_pred_ransac) ** 2)):.3f}")
    print(f"  Inliers:   {inlier_mask.sum()}/{len(y)}")

    # 6. Theil-Sen (median-based, very robust)
    theil = TheilSenRegressor(random_state=42)
    theil.fit(x.reshape(-1, 1), y)
    y_pred_theil = theil.predict(x.reshape(-1, 1))
    results["Theil-Sen"] = theil
    print(f"\nTheil-Sen:")
    print(f"  Intercept: {theil.intercept_:.3f}")
    print(f"  Slope:     {theil.coef_[0]:.3f}")
    print(f"  RMSE:      {np.sqrt(np.mean((y - y_pred_theil) ** 2)):.3f}")

    return results


def plot_comparison(x, y, is_outlier, results):
    """Plot all regression methods"""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot data
    ax.scatter(x[is_outlier == 0], y[is_outlier == 0], c="blue", s=50, alpha=0.6, label="Inliers")
    ax.scatter(
        x[is_outlier == 1],
        y[is_outlier == 1],
        c="red",
        s=100,
        alpha=0.8,
        marker="x",
        linewidths=2,
        label="Outliers",
    )

    # Plot regression lines
    x_plot = np.linspace(x.min(), x.max(), 100)

    colors = ["black", "green", "orange", "purple", "brown", "pink"]
    linestyles = ["-", "--", "-.", ":", "-", "--"]

    # True relationship
    ax.plot(x_plot, 3 + 2 * x_plot, "r--", lw=2, alpha=0.5, label="True (y=3+2x)")

    # OLS
    ols = results["OLS"]
    X_plot = sm.add_constant(x_plot)
    y_ols = ols.predict(X_plot)
    ax.plot(x_plot, y_ols, color=colors[0], linestyle=linestyles[0], lw=2, label="OLS")

    # RLM methods
    rlm_huber = results["RLM (Huber)"]
    y_rlm_huber = rlm_huber.predict(X_plot)
    ax.plot(
        x_plot, y_rlm_huber, color=colors[1], linestyle=linestyles[1], lw=2, label="RLM (Huber)"
    )

    rlm_tukey = results["RLM (Tukey)"]
    y_rlm_tukey = rlm_tukey.predict(X_plot)
    ax.plot(
        x_plot, y_rlm_tukey, color=colors[2], linestyle=linestyles[2], lw=2, label="RLM (Tukey)"
    )

    # sklearn methods
    huber = results["Huber (sklearn)"]
    y_huber = huber.predict(x_plot.reshape(-1, 1))
    ax.plot(
        x_plot, y_huber, color=colors[3], linestyle=linestyles[3], lw=2, label="Huber (sklearn)"
    )

    ransac = results["RANSAC"]
    y_ransac = ransac.predict(x_plot.reshape(-1, 1))
    ax.plot(x_plot, y_ransac, color=colors[4], linestyle=linestyles[4], lw=2, label="RANSAC")

    theil = results["Theil-Sen"]
    y_theil = theil.predict(x_plot.reshape(-1, 1))
    ax.plot(x_plot, y_theil, color=colors[5], linestyle=linestyles[5], lw=2, label="Theil-Sen")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Comparison of Regression Methods with Outliers")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_weights(x, y, is_outlier, results):
    """Plot observation weights for robust methods"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # RLM Huber weights
    ax = axes[0, 0]
    weights = results["RLM (Huber)"].weights
    colors = ["red" if w < 0.5 else "blue" for w in weights]
    ax.scatter(x, weights, c=colors, s=50, alpha=0.6)
    ax.axhline(y=0.5, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel("Weight")
    ax.set_title("RLM (Huber) - Observation Weights")
    ax.grid(True, alpha=0.3)

    # RLM Tukey weights
    ax = axes[0, 1]
    weights = results["RLM (Tukey)"].weights
    colors = ["red" if w < 0.5 else "blue" for w in weights]
    ax.scatter(x, weights, c=colors, s=50, alpha=0.6)
    ax.axhline(y=0.5, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel("Weight")
    ax.set_title("RLM (Tukey) - Observation Weights")
    ax.grid(True, alpha=0.3)

    # RANSAC inlier mask
    ax = axes[1, 0]
    inlier_mask = results["RANSAC"].inlier_mask_
    colors_ransac = ["blue" if m else "red" for m in inlier_mask]
    ax.scatter(x, y, c=colors_ransac, s=50, alpha=0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("RANSAC - Inliers (blue) vs Outliers (red)")
    ax.grid(True, alpha=0.3)

    # True outliers
    ax = axes[1, 1]
    colors_true = ["blue" if o == 0 else "red" for o in is_outlier]
    ax.scatter(x, y, c=colors_true, s=50, alpha=0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("True Labels - Inliers (blue) vs Outliers (red)")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def main():
    """Main execution function"""
    print("Generating data with outliers...")
    x, y, is_outlier = generate_data_with_outliers(n=100, n_outliers=15, seed=42)

    print(f"Dataset: {len(x)} points")
    print(f"  Inliers:  {(is_outlier == 0).sum()}")
    print(f"  Outliers: {(is_outlier == 1).sum()}")
    print(f"True model: y = 3 + 2*x + noise\n")

    # Compare methods
    results = compare_regression_methods(x, y)

    # Plot comparison
    fig1 = plot_comparison(x, y, is_outlier, results)

    # Plot weights
    fig2 = plot_weights(x, y, is_outlier, results)

    # Method recommendations
    print("\n" + "=" * 80)
    print("WHEN TO USE EACH METHOD")
    print("=" * 80)
    print("""
OLS (Ordinary Least Squares):
  - Use when: No outliers present, assumptions met
  - Pros: Simple, interpretable, statistical inference
  - Cons: Very sensitive to outliers

RLM (Huber):
  - Use when: Moderate outliers, want statistical inference
  - Pros: Good balance of robustness and efficiency
  - Cons: Less robust than RANSAC/Theil-Sen

RLM (Tukey):
  - Use when: More aggressive outlier downweighting needed
  - Pros: More robust than Huber
  - Cons: Can be too aggressive

Huber Regressor (sklearn):
  - Use when: Similar to RLM Huber, but using sklearn pipeline
  - Pros: Integrates with sklearn, robust
  - Cons: No p-values/statistical inference

RANSAC:
  - Use when: Many outliers (up to 50%), need very robust fit
  - Pros: Very robust, identifies inlier subset explicitly
  - Cons: Randomized (set seed!), no inference

Theil-Sen:
  - Use when: Small datasets, need median-based robust fit
  - Pros: Very robust, deterministic
  - Cons: Slow for large datasets (O(n²))
    """)

    # Show plots
    print("COMPLETE - Displaying plots")
    plt.show()


if __name__ == "__main__":
    main()

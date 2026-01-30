"""
Complete OLS Regression Workflow with Diagnostics

This example demonstrates:
- Loading and exploring data
- Fitting an OLS regression model
- Comprehensive assumption checking
- Identifying influential points and outliers
- Model refinement
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot
from statsmodels.stats.outliers_influence import variance_inflation_factor, OLSInfluence
from statsmodels.stats.diagnostic import het_breuschpagan, linear_harvey_collier
from statsmodels.stats.stattools import durbin_watson
from scipy.stats import shapiro


def generate_sample_data(n=200, seed=42):
    """Generate sample data for demonstration"""
    np.random.seed(seed)

    # Generate predictors
    x1 = np.random.normal(10, 2, n)
    x2 = np.random.normal(5, 1, n)
    x3 = np.random.normal(20, 3, n)

    # Generate response with some noise
    noise = np.random.normal(0, 2, n)
    y = 3 + 2 * x1 - 1.5 * x2 + 0.5 * x3 + noise

    # Add a few outliers
    outlier_indices = np.random.choice(n, size=5, replace=False)
    y[outlier_indices] += np.random.choice([-15, 15], size=5)

    # Create DataFrame
    df = pd.DataFrame({"y": y, "x1": x1, "x2": x2, "x3": x3})

    return df


def exploratory_analysis(df):
    """Perform exploratory data analysis"""
    print("=" * 80)
    print("EXPLORATORY DATA ANALYSIS")
    print("=" * 80)

    # Summary statistics
    print("\nSummary Statistics:")
    print(df.describe())

    # Correlation matrix
    print("\nCorrelation Matrix:")
    print(df.corr())

    # Visualizations
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Histograms
    df.hist(bins=30, ax=axes.ravel()[:4], edgecolor="black")
    plt.suptitle("Distributions of Variables")
    plt.tight_layout()

    # Pairplot-style scatter matrix
    fig2 = plt.figure(figsize=(12, 10))
    pd.plotting.scatter_matrix(df, figsize=(12, 10), diagonal="hist")
    plt.suptitle("Scatter Matrix")
    plt.tight_layout()

    return fig, fig2


def fit_ols_model(df):
    """Fit OLS regression model"""
    print("\n" + "=" * 80)
    print("FITTING OLS MODEL")
    print("=" * 80)

    # Prepare data
    y = df["y"]
    X = df[["x1", "x2", "x3"]]
    X = sm.add_constant(X)  # Add intercept

    # Fit model
    model = sm.OLS(y, X)
    results = model.fit()

    # Print summary
    print(results.summary())

    return results, X, y


def check_multicollinearity(X):
    """Check for multicollinearity using VIF"""
    print("\n" + "=" * 80)
    print("MULTICOLLINEARITY CHECK (VIF)")
    print("=" * 80)

    vif_data = pd.DataFrame()
    vif_data["Variable"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]

    print("\n", vif_data)
    print("\nInterpretation:")
    print("  VIF > 10: Serious multicollinearity")
    print("  VIF > 5: Moderate multicollinearity")
    print("  VIF < 5: Acceptable")

    return vif_data


def check_assumptions(results):
    """Check regression assumptions with statistical tests"""
    print("\n" + "=" * 80)
    print("ASSUMPTION TESTS")
    print("=" * 80)

    # 1. Linearity (Harvey-Collier test)
    try:
        linearity_test = linear_harvey_collier(results)
        print(f"\n1. Linearity (Harvey-Collier test)")
        print(f"   Test statistic: {linearity_test[0]:.4f}")
        print(f"   P-value: {linearity_test[1]:.4f}")
        interp = "Linear relationship OK" if linearity_test[1] > 0.05 else "Non-linearity detected"
        print(f"   Interpretation: {interp}")
    except:
        print("\n1. Linearity test: Could not be performed (check residual plot manually)")

    # 2. Independence (Durbin-Watson)
    dw = durbin_watson(results.resid)
    print(f"\n2. Independence (Durbin-Watson)")
    print(f"   Statistic: {dw:.4f}")
    print(f"   Interpretation: ", end="")
    if 1.5 < dw < 2.5:
        print("No significant autocorrelation")
    elif dw < 1.5:
        print("Positive autocorrelation detected")
    else:
        print("Negative autocorrelation detected")

    # 3. Normality (Shapiro-Wilk)
    stat, pval = shapiro(results.resid)
    print(f"\n3. Normality of Residuals (Shapiro-Wilk)")
    print(f"   Test statistic: {stat:.4f}")
    print(f"   P-value: {pval:.4f}")
    interp = "Residuals appear normal" if pval > 0.05 else "Non-normality detected"
    print(f"   Interpretation: {interp}")

    # 4. Homoskedasticity (Breusch-Pagan)
    lm_stat, lm_pval, f_stat, f_pval = het_breuschpagan(results.resid, results.model.exog)
    print(f"\n4. Homoskedasticity (Breusch-Pagan)")
    print(f"   LM statistic: {lm_stat:.4f}")
    print(f"   P-value: {lm_pval:.4f}")
    interp = "Constant variance OK" if lm_pval > 0.05 else "Heteroskedasticity detected"
    print(f"   Interpretation: {interp}")

    return {"durbin_watson": dw, "shapiro_pval": pval, "breusch_pagan_pval": lm_pval}


def plot_diagnostics(results):
    """Create comprehensive diagnostic plots"""
    print("\n" + "=" * 80)
    print("DIAGNOSTIC PLOTS")
    print("=" * 80)
    print("Creating diagnostic plots...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Get influence measures
    influence = OLSInfluence(results)

    # 1. Residuals vs Fitted
    axes[0, 0].scatter(results.fittedvalues, results.resid, alpha=0.6)
    axes[0, 0].axhline(y=0, color="r", linestyle="--", linewidth=2)
    axes[0, 0].set_xlabel("Fitted values")
    axes[0, 0].set_ylabel("Residuals")
    axes[0, 0].set_title("Residuals vs Fitted\n(Should show random scatter)")

    # Add lowess smooth
    from statsmodels.nonparametric.smoothers_lowess import lowess

    smoothed = lowess(results.resid, results.fittedvalues, frac=0.3)
    axes[0, 0].plot(smoothed[:, 0], smoothed[:, 1], "g-", lw=2, label="Lowess smooth")
    axes[0, 0].legend()

    # 2. Q-Q Plot
    qqplot(results.resid, line="s", ax=axes[0, 1])
    axes[0, 1].set_title("Normal Q-Q\n(Points should follow line)")

    # 3. Scale-Location
    resid_sqrt = np.sqrt(np.abs(results.resid_pearson))
    axes[1, 0].scatter(results.fittedvalues, resid_sqrt, alpha=0.6)
    axes[1, 0].set_xlabel("Fitted values")
    axes[1, 0].set_ylabel("√|Standardized Residuals|")
    axes[1, 0].set_title("Scale-Location\n(Should show constant spread)")
    smoothed = lowess(resid_sqrt, results.fittedvalues, frac=0.3)
    axes[1, 0].plot(smoothed[:, 0], smoothed[:, 1], "r-", lw=2)

    # 4. Residuals vs Leverage
    leverage = influence.hat_matrix_diag
    cooks_d = influence.cooks_distance[0]
    n = len(results.resid)

    axes[1, 1].scatter(leverage, results.resid_pearson, alpha=0.6)
    axes[1, 1].set_xlabel("Leverage")
    axes[1, 1].set_ylabel("Standardized Residuals")
    axes[1, 1].set_title("Residuals vs Leverage\n(Cook's distance contours)")

    # Highlight high Cook's D points
    high_cooks = cooks_d > 4 / n
    axes[1, 1].scatter(
        leverage[high_cooks],
        results.resid_pearson[high_cooks],
        color="red",
        s=100,
        edgecolor="black",
        linewidth=2,
        label=f"High Cook's D (>{4 / n:.4f})",
    )
    axes[1, 1].legend()

    plt.tight_layout()

    return fig


def identify_influential_points(results, X, y):
    """Identify outliers and influential points"""
    print("\n" + "=" * 80)
    print("INFLUENTIAL POINTS AND OUTLIERS")
    print("=" * 80)

    influence = OLSInfluence(results)
    n, p = X.shape

    # Get influence measures
    cooks_d = influence.cooks_distance[0]
    leverage = influence.hat_matrix_diag
    dffits = influence.dffits[0]
    student_resid = influence.resid_studentized_internal

    # Identify problematic points
    influential = cooks_d > 4 / n
    high_leverage = leverage > 2 * p / n
    outliers = np.abs(student_resid) > 3

    print(f"\nThresholds:")
    print(f"  Cook's D > {4 / n:.4f}")
    print(f"  Leverage > {2 * p / n:.4f}")
    print(f"  |Studentized Residual| > 3")

    print(f"\nSummary:")
    print(f"  Influential points (high Cook's D): {influential.sum()}")
    print(f"  High leverage points: {high_leverage.sum()}")
    print(f"  Outliers in Y-space: {outliers.sum()}")

    # Create detailed table for problematic points
    problematic = influential | high_leverage | outliers

    if problematic.sum() > 0:
        print(f"\nDetailed Information on {problematic.sum()} Problematic Points:")
        influence_df = pd.DataFrame(
            {
                "Index": np.where(problematic)[0],
                "Cooks_D": cooks_d[problematic],
                "Leverage": leverage[problematic],
                "Studentized_Resid": student_resid[problematic],
                "DFFITS": dffits[problematic],
                "y_actual": y[problematic],
                "y_fitted": results.fittedvalues[problematic],
            }
        )
        print(influence_df.to_string())

    return influence, influential


def compare_with_without_influential(df, results, influential):
    """Compare model with and without influential points"""
    print("\n" + "=" * 80)
    print("MODEL COMPARISON: WITH VS WITHOUT INFLUENTIAL POINTS")
    print("=" * 80)

    # Original model metrics
    print("\nOriginal Model (all data):")
    print(f"  R²: {results.rsquared:.4f}")
    print(f"  Adjusted R²: {results.rsquared_adj:.4f}")
    print(f"  RMSE: {np.sqrt(np.mean(results.resid**2)):.4f}")
    print(f"  AIC: {results.aic:.2f}")
    print(f"  BIC: {results.bic:.2f}")

    if influential.sum() > 0:
        # Refit without influential points
        df_clean = df[~influential].copy()
        y_clean = df_clean["y"]
        X_clean = df_clean[["x1", "x2", "x3"]]
        X_clean = sm.add_constant(X_clean)

        model_clean = sm.OLS(y_clean, X_clean)
        results_clean = model_clean.fit()

        print(f"\nModel without {influential.sum()} influential points:")
        print(f"  R²: {results_clean.rsquared:.4f}")
        print(f"  Adjusted R²: {results_clean.rsquared_adj:.4f}")
        print(f"  RMSE: {np.sqrt(np.mean(results_clean.resid**2)):.4f}")
        print(f"  AIC: {results_clean.aic:.2f}")
        print(f"  BIC: {results_clean.bic:.2f}")

        print("\nCoefficient Comparison:")
        coef_comparison = pd.DataFrame(
            {
                "Original": results.params,
                "Without_Influential": results_clean.params,
                "Difference": results.params - results_clean.params,
                "Pct_Change": 100 * (results.params - results_clean.params) / results.params,
            }
        )
        print(coef_comparison.to_string())

        return results_clean
    else:
        print("\nNo influential points to remove.")
        return None


def main():
    """Main execution function"""
    # Generate sample data
    print("Generating sample data...")
    df = generate_sample_data(n=200, seed=42)

    # Exploratory analysis
    fig1, fig2 = exploratory_analysis(df)

    # Fit OLS model
    results, X, y = fit_ols_model(df)

    # Check multicollinearity
    vif_data = check_multicollinearity(X)

    # Check assumptions
    assumption_tests = check_assumptions(results)

    # Create diagnostic plots
    fig_diag = plot_diagnostics(results)

    # Identify influential points
    influence, influential = identify_influential_points(results, X, y)

    # Compare with/without influential points
    results_clean = compare_with_without_influential(df, results, influential)

    # Show all plots
    print("\n" + "=" * 80)
    print("COMPLETE - Displaying all plots")
    print("=" * 80)
    plt.show()


if __name__ == "__main__":
    main()

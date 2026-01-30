"""
Machine Learning Regression Model Comparison

This example demonstrates:
- Comparing multiple ML regression algorithms
- Cross-validation for model selection
- Hyperparameter tuning with GridSearchCV
- Feature importance analysis
- Model evaluation on test set
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, cross_validate, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, HuberRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor,
    ExtraTreesRegressor,
)
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import (
    mean_squared_error,
    mean_absolute_error,
    r2_score,
    mean_absolute_percentage_error,
)


def generate_sample_data(n=1000, n_features=10, seed=42):
    """Generate sample regression data"""
    np.random.seed(seed)

    # Generate features
    X = np.random.randn(n, n_features)

    # True coefficients (only some features are relevant)
    true_coef = np.array([2.5, -1.5, 3.0, 0, 0, 1.2, 0, 0, -0.8, 0.5])

    # Generate response
    y = X @ true_coef + np.random.randn(n) * 0.5

    # Add some nonlinearity
    y += 0.1 * X[:, 0] ** 2 - 0.05 * X[:, 2] ** 2

    # Create feature names
    feature_names = [f"feature_{i + 1}" for i in range(n_features)]

    return X, y, feature_names, true_coef


def compare_models_cv(X, y, cv=5):
    """Compare multiple models using cross-validation"""
    print("=" * 80)
    print("MODEL COMPARISON WITH CROSS-VALIDATION")
    print("=" * 80)

    # Define models to compare
    models = {
        "Linear Regression": LinearRegression(),
        "Ridge (α=1.0)": Ridge(alpha=1.0),
        "Lasso (α=0.1)": Lasso(alpha=0.1),
        "ElasticNet": ElasticNet(alpha=0.1, l1_ratio=0.5),
        "Huber (robust)": HuberRegressor(),
        "KNN (k=5)": KNeighborsRegressor(n_neighbors=5),
        "Random Forest": RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1),
        "Gradient Boosting": GradientBoostingRegressor(n_estimators=100, random_state=42),
        "Extra Trees": ExtraTreesRegressor(n_estimators=100, random_state=42, n_jobs=-1),
    }

    # Store results
    results = []

    print(f"\nRunning {cv}-fold cross-validation for {len(models)} models...\n")

    for name, model in models.items():
        # Create pipeline with scaling
        pipeline = Pipeline([("scaler", StandardScaler()), ("model", model)])

        # Cross-validation with multiple metrics
        cv_results = cross_validate(
            pipeline,
            X,
            y,
            cv=cv,
            scoring=["neg_mean_squared_error", "r2", "neg_mean_absolute_error"],
            return_train_score=True,
            n_jobs=-1,
        )

        # Calculate statistics
        test_rmse = np.sqrt(-cv_results["test_neg_mean_squared_error"])
        train_rmse = np.sqrt(-cv_results["train_neg_mean_squared_error"])
        test_r2 = cv_results["test_r2"]
        test_mae = -cv_results["test_neg_mean_absolute_error"]

        results.append(
            {
                "Model": name,
                "Test_RMSE_mean": test_rmse.mean(),
                "Test_RMSE_std": test_rmse.std(),
                "Train_RMSE_mean": train_rmse.mean(),
                "Test_R2_mean": test_r2.mean(),
                "Test_R2_std": test_r2.std(),
                "Test_MAE_mean": test_mae.mean(),
                "Overfit_Gap": train_rmse.mean() - test_rmse.mean(),
            }
        )

    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("Test_RMSE_mean")

    print("Cross-Validation Results (sorted by RMSE):")
    print("-" * 80)
    print(results_df.to_string(index=False))

    return results_df, models


def tune_best_model(X_train, y_train, cv=5):
    """Hyperparameter tuning for Random Forest"""
    print("\n" + "=" * 80)
    print("HYPERPARAMETER TUNING (Random Forest)")
    print("=" * 80)

    # Define parameter grid
    param_grid = {
        "model__n_estimators": [50, 100, 200],
        "model__max_depth": [5, 10, 15, None],
        "model__min_samples_split": [2, 5, 10],
        "model__min_samples_leaf": [1, 2, 4],
    }

    # Create pipeline
    pipeline = Pipeline(
        [("scaler", StandardScaler()), ("model", RandomForestRegressor(random_state=42, n_jobs=-1))]
    )

    # Grid search
    print(f"\nSearching {np.prod([len(v) for v in param_grid.values()])} combinations...")
    grid_search = GridSearchCV(
        pipeline, param_grid, cv=cv, scoring="neg_mean_squared_error", n_jobs=-1, verbose=1
    )

    grid_search.fit(X_train, y_train)

    print(f"\nBest parameters: {grid_search.best_params_}")
    print(f"Best CV RMSE: {np.sqrt(-grid_search.best_score_):.4f}")

    # Show top 5 parameter combinations
    results_df = pd.DataFrame(grid_search.cv_results_)
    results_df["RMSE"] = np.sqrt(-results_df["mean_test_score"])
    results_df = results_df.sort_values("RMSE")

    print("\nTop 5 parameter combinations:")
    cols = [
        "param_model__n_estimators",
        "param_model__max_depth",
        "param_model__min_samples_split",
        "param_model__min_samples_leaf",
        "RMSE",
    ]
    print(results_df[cols].head().to_string(index=False))

    return grid_search.best_estimator_


def evaluate_on_test_set(model, X_train, X_test, y_train, y_test):
    """Evaluate tuned model on test set"""
    print("\n" + "=" * 80)
    print("TEST SET EVALUATION")
    print("=" * 80)

    # Train on full training set
    model.fit(X_train, y_train)

    # Predictions
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # Calculate metrics
    print("\nTraining Set Performance:")
    print(f"  RMSE: {np.sqrt(mean_squared_error(y_train, y_train_pred)):.4f}")
    print(f"  MAE:  {mean_absolute_error(y_train, y_train_pred):.4f}")
    print(f"  R²:   {r2_score(y_train, y_train_pred):.4f}")

    print("\nTest Set Performance:")
    print(f"  RMSE: {np.sqrt(mean_squared_error(y_test, y_test_pred)):.4f}")
    print(f"  MAE:  {mean_absolute_error(y_test, y_test_pred):.4f}")
    print(f"  R²:   {r2_score(y_test, y_test_pred):.4f}")

    # Calculate residuals
    residuals = y_test - y_test_pred

    print("\nResidual Statistics:")
    print(f"  Mean:   {np.mean(residuals):.4f}")
    print(f"  Std:    {np.std(residuals):.4f}")
    print(f"  Min:    {np.min(residuals):.4f}")
    print(f"  Max:    {np.max(residuals):.4f}")

    return y_test_pred, residuals


def analyze_feature_importance(model, feature_names):
    """Analyze and plot feature importance"""
    print("\n" + "=" * 80)
    print("FEATURE IMPORTANCE ANALYSIS")
    print("=" * 80)

    # Get the actual model from the pipeline
    rf_model = model.named_steps["model"]

    # Get feature importance
    importance = rf_model.feature_importances_

    # Create DataFrame
    importance_df = pd.DataFrame({"Feature": feature_names, "Importance": importance}).sort_values(
        "Importance", ascending=False
    )

    print("\nFeature Importance:")
    print(importance_df.to_string(index=False))

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(range(len(importance_df)), importance_df["Importance"])
    ax.set_yticks(range(len(importance_df)))
    ax.set_yticklabels(importance_df["Feature"])
    ax.set_xlabel("Importance")
    ax.set_title("Feature Importance (Random Forest)")
    ax.invert_yaxis()
    plt.tight_layout()

    return importance_df, fig


def plot_results(y_test, y_pred, residuals):
    """Create visualization of results"""
    print("\n" + "=" * 80)
    print("CREATING VISUALIZATIONS")
    print("=" * 80)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Actual vs Predicted
    axes[0, 0].scatter(y_test, y_pred, alpha=0.5)
    min_val = min(y_test.min(), y_pred.min())
    max_val = max(y_test.max(), y_pred.max())
    axes[0, 0].plot([min_val, max_val], [min_val, max_val], "r--", lw=2)
    axes[0, 0].set_xlabel("Actual")
    axes[0, 0].set_ylabel("Predicted")
    axes[0, 0].set_title("Actual vs Predicted")
    axes[0, 0].grid(True, alpha=0.3)

    # 2. Residuals vs Predicted
    axes[0, 1].scatter(y_pred, residuals, alpha=0.5)
    axes[0, 1].axhline(y=0, color="r", linestyle="--", lw=2)
    axes[0, 1].set_xlabel("Predicted")
    axes[0, 1].set_ylabel("Residuals")
    axes[0, 1].set_title("Residuals vs Predicted")
    axes[0, 1].grid(True, alpha=0.3)

    # 3. Residual Distribution
    axes[1, 0].hist(residuals, bins=50, edgecolor="black", alpha=0.7)
    axes[1, 0].axvline(x=0, color="r", linestyle="--", lw=2)
    axes[1, 0].set_xlabel("Residuals")
    axes[1, 0].set_ylabel("Frequency")
    axes[1, 0].set_title("Residual Distribution")

    # 4. Q-Q Plot
    from scipy import stats

    stats.probplot(residuals, dist="norm", plot=axes[1, 1])
    axes[1, 1].set_title("Q-Q Plot")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    return fig


def plot_cv_comparison(results_df):
    """Plot cross-validation comparison"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # RMSE comparison
    ax = axes[0]
    y_pos = np.arange(len(results_df))
    ax.barh(
        y_pos, results_df["Test_RMSE_mean"], xerr=results_df["Test_RMSE_std"], alpha=0.7, capsize=5
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df["Model"])
    ax.invert_yaxis()
    ax.set_xlabel("Test RMSE (lower is better)")
    ax.set_title("Cross-Validation RMSE Comparison")
    ax.grid(True, alpha=0.3, axis="x")

    # R² comparison
    ax = axes[1]
    ax.barh(
        y_pos,
        results_df["Test_R2_mean"],
        xerr=results_df["Test_R2_std"],
        alpha=0.7,
        capsize=5,
        color="green",
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df["Model"])
    ax.invert_yaxis()
    ax.set_xlabel("Test R² (higher is better)")
    ax.set_title("Cross-Validation R² Comparison")
    ax.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()

    return fig


def main():
    """Main execution function"""
    print("Generating sample data...")
    X, y, feature_names, true_coef = generate_sample_data(n=1000, n_features=10)

    print(f"Dataset: {X.shape[0]} samples, {X.shape[1]} features")

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    print(f"Training set: {X_train.shape[0]} samples")
    print(f"Test set: {X_test.shape[0]} samples\n")

    # Compare models with cross-validation
    results_df, models = compare_models_cv(X_train, y_train, cv=5)

    # Plot CV comparison
    fig_cv = plot_cv_comparison(results_df)

    # Tune best model (Random Forest)
    best_model = tune_best_model(X_train, y_train, cv=5)

    # Evaluate on test set
    y_pred, residuals = evaluate_on_test_set(best_model, X_train, X_test, y_train, y_test)

    # Feature importance
    importance_df, fig_importance = analyze_feature_importance(best_model, feature_names)

    # Plot results
    fig_results = plot_results(y_test, y_pred, residuals)

    # Show comparison of true vs estimated importance
    print("\n" + "=" * 80)
    print("TRUE VS ESTIMATED FEATURE IMPORTANCE")
    print("=" * 80)

    comparison = pd.DataFrame(
        {
            "Feature": feature_names,
            "True_Coefficient": true_coef,
            "Estimated_Importance": importance_df.set_index("Feature")
            .loc[feature_names, "Importance"]
            .values,
        }
    )
    print("\n", comparison.to_string(index=False))

    # Display all plots
    print("\n" + "=" * 80)
    print("COMPLETE - Displaying all plots")
    print("=" * 80)
    plt.show()


if __name__ == "__main__":
    main()

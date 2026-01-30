"""
Complete Experimental Workflow with Claude-Light

This script demonstrates a comprehensive experimental workflow including:
1. Statistics and reproducibility
2. Linear regression
3. Multivariate modeling
4. Optimization
5. Design of experiments
6. Machine learning
"""

import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import minimize
from scipy.stats import qmc
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score

# API endpoint
API_URL = "https://claude-light.cheme.cmu.edu/api"


def get_measurement(R, G, B, timeout=10):
    """Get measurement from Claude-Light."""
    resp = requests.get(API_URL, params={"R": R, "G": G, "B": B}, timeout=timeout)
    return resp.json()


def experiment_1_statistics():
    """Experiment 1: Assess measurement reproducibility."""
    print("=" * 60)
    print("EXPERIMENT 1: Statistics and Reproducibility")
    print("=" * 60)

    # Fixed conditions
    R, G, B = 0.5, 0.5, 0.5
    n_repeats = 30

    # Collect repeated measurements
    measurements_515 = []
    measurements_630 = []

    print(f"Taking {n_repeats} measurements at R={R}, G={G}, B={B}...")

    for i in range(n_repeats):
        data = get_measurement(R, G, B)
        measurements_515.append(data["out"]["515nm"])
        measurements_630.append(data["out"]["630nm"])
        if (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{n_repeats} measurements")

    # Statistical analysis
    measurements_515 = np.array(measurements_515)
    measurements_630 = np.array(measurements_630)

    print("\nResults for 515nm channel:")
    print(f"  Mean: {np.mean(measurements_515):.2f}")
    print(f"  Std Dev: {np.std(measurements_515):.2f}")
    print(f"  CV: {100 * np.std(measurements_515) / np.mean(measurements_515):.2f}%")
    print(f"  Median: {np.median(measurements_515):.2f}")
    print(
        f"  95% CI: [{np.percentile(measurements_515, 2.5):.2f}, "
        f"{np.percentile(measurements_515, 97.5):.2f}]"
    )

    # Plot histogram
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.hist(measurements_515, bins=15, edgecolor="black")
    plt.xlabel("Intensity at 515nm")
    plt.ylabel("Count")
    plt.title("Measurement Distribution (515nm)")

    plt.subplot(1, 2, 2)
    plt.hist(measurements_630, bins=15, edgecolor="black", color="red")
    plt.xlabel("Intensity at 630nm")
    plt.ylabel("Count")
    plt.title("Measurement Distribution (630nm)")
    plt.tight_layout()
    plt.savefig("experiment1_statistics.png", dpi=150)
    print("\nSaved histogram to 'experiment1_statistics.png'")
    plt.close()

    return measurements_515, measurements_630


def experiment_2_linear_regression():
    """Experiment 2: Single-variable linear regression."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 2: Linear Regression (Red LED)")
    print("=" * 60)

    # Vary R from 0 to 1
    R_values = np.linspace(0, 1, 11)
    outputs_630 = []

    print("Measuring red LED response at 630nm...")

    for R in R_values:
        data = get_measurement(R, 0, 0)
        outputs_630.append(data["out"]["630nm"])
        print(f"  R = {R:.2f}: Output = {outputs_630[-1]:.2f}")

    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(R_values, outputs_630)

    print(f"\nLinear Fit Results:")
    print(f"  Slope: {slope:.2f}")
    print(f"  Intercept: {intercept:.2f}")
    print(f"  R²: {r_value**2:.4f}")
    print(f"  Std Error: {std_err:.2f}")

    # Prediction
    target_output = 25000
    predicted_R = (target_output - intercept) / slope
    print(f"\nPrediction:")
    print(f"  To achieve output of {target_output}, R should be {predicted_R:.3f}")

    # Validate prediction
    print(f"\nValidating prediction...")
    validation_data = get_measurement(predicted_R, 0, 0)
    actual_output = validation_data["out"]["630nm"]
    error = abs(actual_output - target_output)
    print(f"  Predicted R: {predicted_R:.3f}")
    print(f"  Target output: {target_output}")
    print(f"  Actual output: {actual_output:.2f}")
    print(f"  Error: {error:.2f} ({100 * error / target_output:.2f}%)")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(R_values, outputs_630, "o", label="Measured", markersize=8)
    plt.plot(R_values, slope * R_values + intercept, "r-", label="Linear Fit")
    plt.plot(predicted_R, target_output, "s", markersize=10, label="Prediction", color="green")
    plt.xlabel("R Input")
    plt.ylabel("Output at 630nm")
    plt.title(f"Red LED Response (R² = {r_value**2:.4f})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("experiment2_regression.png", dpi=150)
    print("\nSaved plot to 'experiment2_regression.png'")
    plt.close()

    return R_values, outputs_630, slope, intercept


def experiment_3_multivariate():
    """Experiment 3: Multivariate regression."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 3: Multivariate Regression")
    print("=" * 60)

    # Latin hypercube sampling
    sampler = qmc.LatinHypercube(d=3)
    n_samples = 25
    samples = sampler.random(n=n_samples)
    samples = qmc.scale(samples, [0, 0, 0], [1, 1, 1])

    print(f"Collecting {n_samples} measurements using Latin Hypercube Sampling...")

    X = []  # Inputs (R, G, B)
    y_515 = []  # Output at 515nm
    y_630 = []  # Output at 630nm

    for i, (R, G, B) in enumerate(samples):
        data = get_measurement(R, G, B)
        X.append([R, G, B])
        y_515.append(data["out"]["515nm"])
        y_630.append(data["out"]["630nm"])

        if (i + 1) % 5 == 0:
            print(f"  Completed {i + 1}/{n_samples} measurements")

    X = np.array(X)
    y_515 = np.array(y_515)
    y_630 = np.array(y_630)

    # Fit model for 515nm
    model_515 = LinearRegression()
    model_515.fit(X, y_515)

    print(f"\nModel for 515nm:")
    print(f"  R coefficient: {model_515.coef_[0]:.2f}")
    print(f"  G coefficient: {model_515.coef_[1]:.2f}")
    print(f"  B coefficient: {model_515.coef_[2]:.2f}")
    print(f"  Intercept: {model_515.intercept_:.2f}")
    print(f"  R² score: {model_515.score(X, y_515):.4f}")

    # Fit model for 630nm
    model_630 = LinearRegression()
    model_630.fit(X, y_630)

    print(f"\nModel for 630nm:")
    print(f"  R coefficient: {model_630.coef_[0]:.2f}")
    print(f"  G coefficient: {model_630.coef_[1]:.2f}")
    print(f"  B coefficient: {model_630.coef_[2]:.2f}")
    print(f"  Intercept: {model_630.intercept_:.2f}")
    print(f"  R² score: {model_630.score(X, y_630):.4f}")

    # Save data
    df = pd.DataFrame(X, columns=["R", "G", "B"])
    df["output_515nm"] = y_515
    df["output_630nm"] = y_630
    df.to_csv("experiment3_multivariate_data.csv", index=False)
    print("\nSaved data to 'experiment3_multivariate_data.csv'")

    return X, y_515, y_630, model_515, model_630


def experiment_4_optimization():
    """Experiment 4: Optimization to find target output."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 4: Optimization")
    print("=" * 60)

    # Target: 30000 at 515nm and 20000 at 630nm
    target_515 = 30000
    target_630 = 20000

    print(f"Optimizing to achieve:")
    print(f"  515nm: {target_515}")
    print(f"  630nm: {target_630}")

    call_count = [0]  # Mutable to track in nested function

    def objective(inputs):
        R, G, B = np.clip(inputs, 0, 1)
        data = get_measurement(R, G, B)

        actual_515 = data["out"]["515nm"]
        actual_630 = data["out"]["630nm"]

        error = (actual_515 - target_515) ** 2 + (actual_630 - target_630) ** 2

        call_count[0] += 1
        if call_count[0] % 5 == 0:
            print(
                f"  Iteration {call_count[0]}: R={R:.3f}, G={G:.3f}, B={B:.3f}, Error={error:.0f}"
            )

        return error

    # Optimize
    result = minimize(
        objective,
        [0.5, 0.5, 0.5],
        bounds=[(0, 1), (0, 1), (0, 1)],
        method="Nelder-Mead",
        options={"maxiter": 50},
    )

    print(f"\nOptimization completed in {call_count[0]} function calls")
    print(f"Optimal RGB: R={result.x[0]:.3f}, G={result.x[1]:.3f}, B={result.x[2]:.3f}")

    # Validate
    print("\nValidating optimal solution...")
    validation = get_measurement(*result.x)
    print(f"  Target 515nm: {target_515}, Actual: {validation['out']['515nm']:.2f}")
    print(f"  Target 630nm: {target_630}, Actual: {validation['out']['630nm']:.2f}")

    return result.x


def experiment_5_machine_learning():
    """Experiment 5: Machine learning models."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 5: Machine Learning")
    print("=" * 60)

    # Use data from experiment 3
    print("Reusing data from Experiment 3 for ML models...")

    # Generate more data for ML
    print("Collecting additional data for ML training...")
    n_additional = 25

    X = []
    y = []

    for i in range(n_additional):
        R, G, B = np.random.random(3)
        data = get_measurement(R, G, B)
        X.append([R, G, B])
        y.append(data["out"]["515nm"])

        if (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{n_additional} measurements")

    X = np.array(X)
    y = np.array(y)

    # Test different models
    models = {
        "Linear Regression": LinearRegression(),
        "Ridge Regression": Ridge(alpha=1.0),
        "Random Forest": RandomForestRegressor(n_estimators=50, random_state=42),
    }

    print("\nModel Comparison (5-fold cross-validation):")
    for name, model in models.items():
        scores = cross_val_score(model, X, y, cv=5, scoring="r2")
        print(f"  {name:20s}: R² = {scores.mean():.4f} (+/- {scores.std():.4f})")

    # Train final model
    best_model = RandomForestRegressor(n_estimators=100, random_state=42)
    best_model.fit(X, y)

    # Feature importance
    print(f"\nFeature Importance (Random Forest):")
    print(f"  R: {best_model.feature_importances_[0]:.3f}")
    print(f"  G: {best_model.feature_importances_[1]:.3f}")
    print(f"  B: {best_model.feature_importances_[2]:.3f}")

    return best_model


if __name__ == "__main__":
    print("Claude-Light Complete Experimental Workflow")
    print("=" * 60)
    print("This will run 5 experiments:")
    print("  1. Statistics and Reproducibility")
    print("  2. Linear Regression")
    print("  3. Multivariate Modeling")
    print("  4. Optimization")
    print("  5. Machine Learning")
    print()
    input("Press Enter to start...")

    # Run all experiments
    try:
        # Experiment 1
        meas_515, meas_630 = experiment_1_statistics()

        # Experiment 2
        R_vals, outputs, slope, intercept = experiment_2_linear_regression()

        # Experiment 3
        X, y_515, y_630, model_515, model_630 = experiment_3_multivariate()

        # Experiment 4
        optimal_rgb = experiment_4_optimization()

        # Experiment 5
        ml_model = experiment_5_machine_learning()

        print("\n" + "=" * 60)
        print("ALL EXPERIMENTS COMPLETED SUCCESSFULLY")
        print("=" * 60)
        print("\nGenerated files:")
        print("  - experiment1_statistics.png")
        print("  - experiment2_regression.png")
        print("  - experiment3_multivariate_data.csv")

    except Exception as e:
        print(f"\nError occurred: {e}")
        print("Check your internet connection and try again.")

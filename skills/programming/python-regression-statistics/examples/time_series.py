"""
Time Series Regression and Forecasting with ARIMA

Demonstrates:
- ARIMA model fitting
- Forecasting with confidence intervals
- Model diagnostics
- Seasonal decomposition
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.stats.diagnostic import acorr_ljungbox


def generate_time_series(n=200, seed=42):
    """Generate sample time series data"""
    np.random.seed(seed)

    # Time index
    dates = pd.date_range("2020-01-01", periods=n, freq="D")

    # Trend + seasonal + noise
    t = np.arange(n)
    trend = 0.5 * t
    seasonal = 10 * np.sin(2 * np.pi * t / 30)  # 30-day cycle
    noise = np.random.normal(0, 3, n)

    y = 50 + trend + seasonal + noise

    # Create series
    ts = pd.Series(y, index=dates)

    return ts


def plot_time_series(ts):
    """Plot time series and decomposition"""
    print("=" * 80)
    print("TIME SERIES VISUALIZATION")
    print("=" * 80)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Original series
    axes[0].plot(ts)
    axes[0].set_title("Original Time Series")
    axes[0].set_ylabel("Value")
    axes[0].grid(True, alpha=0.3)

    # Seasonal decomposition
    decomposition = seasonal_decompose(ts, model="additive", period=30)

    # Plot components
    fig2, axes2 = plt.subplots(4, 1, figsize=(12, 10))
    decomposition.observed.plot(ax=axes2[0])
    axes2[0].set_ylabel("Observed")
    axes2[0].grid(True, alpha=0.3)

    decomposition.trend.plot(ax=axes2[1])
    axes2[1].set_ylabel("Trend")
    axes2[1].grid(True, alpha=0.3)

    decomposition.seasonal.plot(ax=axes2[2])
    axes2[2].set_ylabel("Seasonal")
    axes2[2].grid(True, alpha=0.3)

    decomposition.resid.plot(ax=axes2[3])
    axes2[3].set_ylabel("Residual")
    axes2[3].set_xlabel("Date")
    axes2[3].grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, fig2


def plot_acf_pacf(ts):
    """Plot ACF and PACF"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    plot_acf(ts, lags=40, ax=axes[0])
    axes[0].set_title("Autocorrelation Function (ACF)")

    plot_pacf(ts, lags=40, ax=axes[1])
    axes[1].set_title("Partial Autocorrelation Function (PACF)")

    plt.tight_layout()
    return fig


def fit_arima_model(ts, order=(1, 1, 1)):
    """Fit ARIMA model"""
    print("\n" + "=" * 80)
    print(f"FITTING ARIMA{order} MODEL")
    print("=" * 80)

    # Split into train/test
    train_size = int(len(ts) * 0.8)
    train = ts[:train_size]
    test = ts[train_size:]

    print(f"\nTrain size: {len(train)}")
    print(f"Test size:  {len(test)}")

    # Fit model
    model = ARIMA(train, order=order)
    results = model.fit()

    print("\nModel Summary:")
    print(results.summary())

    return results, train, test


def forecast_and_evaluate(results, train, test, steps=None):
    """Generate forecasts and evaluate"""
    print("\n" + "=" * 80)
    print("FORECASTING AND EVALUATION")
    print("=" * 80)

    if steps is None:
        steps = len(test)

    # Forecast
    forecast_obj = results.get_forecast(steps=steps)
    forecast_df = forecast_obj.summary_frame(alpha=0.05)

    # Calculate metrics on test set
    actual = test[:steps]
    predicted = forecast_df["mean"]

    mae = np.mean(np.abs(actual - predicted))
    rmse = np.sqrt(np.mean((actual - predicted) ** 2))
    mape = np.mean(np.abs((actual - predicted) / actual)) * 100

    print(f"\nForecast Performance:")
    print(f"  MAE:  {mae:.3f}")
    print(f"  RMSE: {rmse:.3f}")
    print(f"  MAPE: {mape:.2f}%")

    return forecast_df


def plot_forecast(train, test, forecast_df):
    """Plot forecast with confidence intervals"""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot train and test
    ax.plot(train.index, train.values, label="Train", color="blue")
    ax.plot(test.index, test.values, label="Test (Actual)", color="green")

    # Plot forecast
    forecast_index = test.index[: len(forecast_df)]
    ax.plot(forecast_index, forecast_df["mean"], label="Forecast", color="red", linestyle="--")

    # Plot confidence interval
    ax.fill_between(
        forecast_index,
        forecast_df["mean_ci_lower"],
        forecast_df["mean_ci_upper"],
        alpha=0.3,
        color="red",
        label="95% CI",
    )

    ax.set_xlabel("Date")
    ax.set_ylabel("Value")
    ax.set_title("ARIMA Forecast with Confidence Intervals")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def check_residuals(results):
    """Check model residuals"""
    print("\n" + "=" * 80)
    print("RESIDUAL DIAGNOSTICS")
    print("=" * 80)

    residuals = results.resid

    # Ljung-Box test for autocorrelation
    lb_test = acorr_ljungbox(residuals, lags=10)
    print("\nLjung-Box Test (residual autocorrelation):")
    print("  Lags with p < 0.05 suggest remaining autocorrelation")
    print(lb_test)

    # Plot residuals
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Residuals over time
    axes[0, 0].plot(residuals)
    axes[0, 0].axhline(y=0, color="r", linestyle="--")
    axes[0, 0].set_title("Residuals over Time")
    axes[0, 0].set_ylabel("Residual")
    axes[0, 0].grid(True, alpha=0.3)

    # Histogram
    axes[0, 1].hist(residuals, bins=30, edgecolor="black")
    axes[0, 1].set_title("Residual Distribution")
    axes[0, 1].set_xlabel("Residual")
    axes[0, 1].set_ylabel("Frequency")

    # Q-Q plot
    from scipy import stats

    stats.probplot(residuals, dist="norm", plot=axes[1, 0])
    axes[1, 0].set_title("Q-Q Plot")

    # ACF of residuals
    plot_acf(residuals, lags=20, ax=axes[1, 1])
    axes[1, 1].set_title("ACF of Residuals")

    plt.tight_layout()
    return fig


def main():
    """Main execution function"""
    print("Generating time series data...")
    ts = generate_time_series(n=200, seed=42)

    # Visualize
    fig1, fig2 = plot_time_series(ts)

    # ACF/PACF
    fig3 = plot_acf_pacf(ts)

    # Fit ARIMA model
    results, train, test = fit_arima_model(ts, order=(1, 1, 1))

    # Forecast
    forecast_df = forecast_and_evaluate(results, train, test)

    # Plot forecast
    fig4 = plot_forecast(train, test, forecast_df)

    # Check residuals
    fig5 = check_residuals(results)

    print("\n" + "=" * 80)
    print("COMPLETE - Displaying plots")
    print("=" * 80)
    plt.show()


if __name__ == "__main__":
    main()

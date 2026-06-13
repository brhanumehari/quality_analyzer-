/*
 * Quality Engineering Analyzer - C++ Implementation for Termux
 * Fixed version - all declarations in correct order
 * 
 * Compile: g++ -std=c++17 -O3 -o quality_analyzer quality_analyzer.cpp -lm
 * Run:     ./quality_analyzer
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include <numeric>
#include <functional>
#include <cstdlib>
#include <ctime>
#include <random>

// ============================================================
// MATH UTILITIES
// ============================================================

const double PI = 3.14159265358979323846;
const double E = 2.71828182845904523536;

double gamma_ln(double x) {
    if (x < 0.5) {
        return log(PI / fabs(sin(PI * x))) - gamma_ln(1.0 - x);
    }
    x -= 1.0;
    double g = 7.0;
    double c[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    double a = c[0];
    for (int i = 1; i < 9; i++) a += c[i] / (x + i);
    double t = x + g + 0.5;
    return 0.5 * log(2 * PI) + (x + 0.5) * log(t) - t + log(a);
}

double gamma_fn(double x) {
    return exp(gamma_ln(x));
}

// Continued fraction for incomplete beta (Lentz's method)
double betainc_cf(double a, double b, double x) {
    const int MAX_ITER = 200;
    const double EPS = 1e-14;
    
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    double h = d;
    
    for (int m = 1; m <= MAX_ITER; m++) {
        int m2 = 2 * m;
        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    return h;
}

// Regularized incomplete beta function
double betainc(double a, double b, double x) {
    if (x < 0 || x > 1) return 0;
    if (x == 0 || x == 1) return x;
    
    double bt = exp(gamma_ln(a + b) - gamma_ln(a) - gamma_ln(b) 
                    + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betainc_cf(a, b, x) / a;
    } else {
        return 1.0 - bt * betainc_cf(b, a, 1.0 - x) / b;
    }
}

// Regularized lower incomplete gamma function P(a,x)
double gamma_reg(double a, double x) {
    if (x <= 0 || a <= 0) return 0;
    
    const int MAX_ITER = 200;
    const double EPS = 1e-14;
    
    double sum = 1.0 / a;
    double term = sum;
    for (int n = 1; n <= MAX_ITER; n++) {
        term *= x / (a + n);
        sum += term;
        if (fabs(term) < fabs(sum) * EPS) break;
    }
    return sum * exp(-x + a * log(x) - gamma_ln(a));
}

// Error function
double erf_custom(double x) {
    double sign = (x >= 0) ? 1.0 : -1.0;
    x = fabs(x);
    double t = 1.0 / (1.0 + 0.3275911 * x);
    double y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t 
                    - 0.284496736) * t + 0.254829592) * t * exp(-x * x);
    return sign * y;
}

double norm_cdf(double x) {
    return 0.5 * (1.0 + erf_custom(x / sqrt(2.0)));
}

double norm_ppf(double p) {
    if (p < 0 || p > 1) return 0;
    if (p < 1e-15) return -8.0;
    if (p > 1 - 1e-15) return 8.0;
    
    // Rational approximation
    double q, r;
    if (p < 0.5) {
        q = sqrt(-2.0 * log(p));
        r = (((((2.3212128 * q + 4.850141) * q + 4.850141) * q + 2.3212128) * q + 1.0)
            / (((0.2316419 * q + 2.3132636) * q + 10.0) * q + 1.0));
        return -r;
    } else {
        return -norm_ppf(1.0 - p);
    }
}

double chi2_cdf(double x, int df) {
    if (x <= 0) return 0;
    return gamma_reg(df / 2.0, x / 2.0);
}

double chi2_ppf(double p, int df) {
    if (p <= 0) return 0;
    if (p >= 1) return 1e10;
    
    double x = df + sqrt(2.0 * df) * norm_ppf(p);
    if (x < 0) x = 0.1;
    
    for (int i = 0; i < 50; i++) {
        double fx = chi2_cdf(x, df) - p;
        if (fabs(fx) < 1e-12) break;
        // Derivative = chi2 pdf
        double dfx = exp((df/2.0 - 1.0) * log(x) - x/2.0 - gamma_ln(df/2.0) - log(2.0));
        if (dfx == 0) break;
        x = x - fx / dfx;
        if (x < 0) x = 0.01;
    }
    return x;
}

double t_cdf(double t, int df) {
    double x = df / (df + t * t);
    if (t >= 0) return 1.0 - 0.5 * betainc(df / 2.0, 0.5, x);
    else return 0.5 * betainc(df / 2.0, 0.5, x);
}

double t_ppf(double p, int df) {
    double x = norm_ppf(p);
    // Hill's approximation for t distribution
    double g1 = (x*x*x + x) / 4.0;
    double g2 = (5*x*x*x*x*x + 16*x*x*x + 3*x) / 96.0;
    double g3 = (3*x*x*x*x*x*x*x + 19*x*x*x*x*x + 17*x*x*x - 15*x) / 384.0;
    return x + g1/df + g2/(df*df) + g3/(df*df*df);
}

// ============================================================
// PDF FUNCTIONS
// ============================================================

double normal_pdf(double x, double mean, double stddev) {
    double z = (x - mean) / stddev;
    return exp(-0.5 * z * z) / (stddev * sqrt(2.0 * PI));
}

double lognormal_pdf(double x, double shape, double loc, double scale) {
    if (x <= loc) return 0;
    double y = (log(x - loc) - scale) / shape;
    return exp(-0.5 * y * y) / ((x - loc) * shape * sqrt(2.0 * PI));
}

double expon_pdf(double x, double loc, double scale) {
    if (x < loc) return 0;
    return exp(-(x - loc) / scale) / scale;
}

double gamma_pdf(double x, double shape, double loc, double scale) {
    if (x < loc) return 0;
    double y = (x - loc) / scale;
    if (y < 0) return 0;
    return pow(y, shape - 1) * exp(-y) / (scale * gamma_fn(shape));
}

double weibull_pdf(double x, double shape, double loc, double scale) {
    if (x < loc) return 0;
    double y = (x - loc) / scale;
    if (y < 0) return 0;
    return (shape / scale) * pow(y, shape - 1) * exp(-pow(y, shape));
}

// ============================================================
// CDF FUNCTIONS FOR DISTRIBUTIONS
// ============================================================

std::function<double(double)> make_normal_cdf(double mean, double stddev) {
    return [mean, stddev](double x) {
        return norm_cdf((x - mean) / stddev);
    };
}

std::function<double(double)> make_lognormal_cdf(double shape, double loc, double scale) {
    return [shape, loc, scale](double x) {
        if (x <= loc) return 0.0;
        return norm_cdf((log(x - loc) - scale) / shape);
    };
}

std::function<double(double)> make_expon_cdf(double loc, double scale) {
    return [loc, scale](double x) {
        if (x < loc) return 0.0;
        return 1.0 - exp(-(x - loc) / scale);
    };
}

std::function<double(double)> make_gamma_cdf(double shape, double loc, double scale) {
    return [shape, loc, scale](double x) {
        if (x < loc) return 0.0;
        return gamma_reg(shape, (x - loc) / scale);
    };
}

std::function<double(double)> make_weibull_cdf(double shape, double loc, double scale) {
    return [shape, loc, scale](double x) {
        if (x < loc) return 0.0;
        return 1.0 - exp(-pow((x - loc) / scale, shape));
    };
}

// ============================================================
// KOLMOGOROV-SMIRNOV TEST
// ============================================================

double ks_test_statistic(const std::vector<double>& data, 
                          std::function<double(double)> cdf) {
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    int n = sorted.size();
    
    double d = 0.0;
    for (int i = 0; i < n; i++) {
        double ecdf = (i + 1.0) / n;
        double cdf_val = cdf(sorted[i]);
        double diff1 = fabs(ecdf - cdf_val);
        double diff2 = fabs((double)i / n - cdf_val);
        d = std::max({d, diff1, diff2});
    }
    return d;
}

double ks_p_value(double d, int n) {
    double sqrt_n = sqrt((double)n);
    double lambda = (sqrt_n + 0.12 + 0.11 / sqrt_n) * d;
    
    // Kolmogorov's asymptotic distribution
    double p = 0.0;
    for (int k = 1; k <= 100; k++) {
        double term = exp(-(2*k-1)*(2*k-1) * PI * PI / (8 * lambda * lambda));
        if (k % 2 == 0) p -= term;
        else p += term;
    }
    p *= sqrt(2.0 * PI) / lambda;
    p = std::max(0.0, std::min(1.0, p));
    return p;
}

// ============================================================
// SHAPIRO-WILK TEST (simplified)
// ============================================================

double shapiro_wilk_p(const std::vector<double>& data) {
    int n = data.size();
    if (n < 3 || n > 5000) return -1;
    
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;
    double s2 = 0;
    for (auto v : data) s2 += (v - mean) * (v - mean);
    
    if (s2 == 0) return 1.0;
    
    // Simplified W statistic
    double b = 0;
    int m = n / 2;
    
    for (int i = 0; i < m; i++) {
        double u = norm_ppf((i + 1.0) / (n + 1.0));
        double a = u / sqrt((double)(n + 1) * (n + 2));
        b += a * (sorted[n - 1 - i] - sorted[i]);
    }
    
    double w = b * b / s2;
    
    // Approximate p-value
    double mu = 0.5;
    double sigma = 1.0 / sqrt(3.0 * (n + 7.0));
    double z = (w - mu) / sigma;
    return norm_cdf(z);
}

// ============================================================
// ANDERSON-DARLING TEST
// ============================================================

double anderson_darling_stat(const std::vector<double>& data) {
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    int n = sorted.size();
    
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;
    double stddev = sqrt(std::accumulate(data.begin(), data.end(), 0.0,
        [mean](double s, double v) { return s + (v - mean) * (v - mean); }) / (n - 1));
    
    double a2 = 0;
    for (int i = 0; i < n; i++) {
        double cdf_val = norm_cdf((sorted[i] - mean) / stddev);
        if (cdf_val <= 0 || cdf_val >= 1) continue;
        a2 += (2.0 * (i + 1) - 1.0) * (log(cdf_val) + log(1.0 - cdf_val));
    }
    a2 = -n - a2 / n;
    
    // Modified statistic for normal distribution
    return a2 * (1.0 + 0.75 / n + 2.25 / (n * n));
}

double anderson_darling_p(double a2_mod) {
    if (a2_mod < 0.576) return 0.15;
    if (a2_mod < 0.656) return 0.10;
    if (a2_mod < 0.787) return 0.05;
    if (a2_mod < 0.918) return 0.025;
    if (a2_mod < 1.092) return 0.01;
    return 0.005;
}

// ============================================================
// DISTRIBUTION FITTING FUNCTIONS
// ============================================================

struct DistFitResult {
    std::string name;
    double p_value;
    double ks_stat;
    std::vector<double> params;
    bool valid;
};

std::pair<double, double> fit_normal(const std::vector<double>& data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double var = std::accumulate(data.begin(), data.end(), 0.0,
        [mean](double s, double v) { return s + (v - mean) * (v - mean); }) / data.size();
    return {mean, sqrt(var)};
}

std::tuple<double, double, double> fit_lognormal(const std::vector<double>& data) {
    double min_val = *std::min_element(data.begin(), data.end());
    double loc = std::max(0.0, min_val - 0.001);
    std::vector<double> log_data;
    for (auto v : data) {
        if (v > loc) log_data.push_back(log(v - loc));
    }
    if (log_data.size() < 2) return {0.1, 0, 0};
    double mean_log = std::accumulate(log_data.begin(), log_data.end(), 0.0) / log_data.size();
    double var_log = std::accumulate(log_data.begin(), log_data.end(), 0.0,
        [mean_log](double s, double v) { return s + (v - mean_log) * (v - mean_log); }) / log_data.size();
    return {sqrt(var_log), loc, mean_log};
}

std::pair<double, double> fit_exponential(const std::vector<double>& data) {
    double min_val = *std::min_element(data.begin(), data.end());
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    return {min_val, mean - min_val};
}

std::tuple<double, double, double> fit_gamma(const std::vector<double>& data) {
    double min_val = *std::min_element(data.begin(), data.end());
    double loc = std::max(0.0, min_val - 0.001);
    std::vector<double> shifted;
    for (auto v : data) shifted.push_back(v - loc);
    
    double mean = std::accumulate(shifted.begin(), shifted.end(), 0.0) / shifted.size();
    double var = std::accumulate(shifted.begin(), shifted.end(), 0.0,
        [mean](double s, double v) { return s + (v - mean) * (v - mean); }) / shifted.size();
    
    if (mean <= 0 || var <= 0) return {1.0, loc, mean};
    
    double shape = mean * mean / var;
    double scale = var / mean;
    return {shape, loc, scale};
}

std::tuple<double, double, double> fit_weibull(const std::vector<double>& data) {
    double min_val = *std::min_element(data.begin(), data.end());
    double loc = std::max(0.0, min_val - 0.001);
    std::vector<double> shifted;
    for (auto v : data) shifted.push_back(v - loc);
    
    int n = shifted.size();
    double mean = std::accumulate(shifted.begin(), shifted.end(), 0.0) / n;
    
    double cv = sqrt(std::accumulate(shifted.begin(), shifted.end(), 0.0,
        [mean](double s, double v) { return s + (v - mean) * (v - mean); }) / (n - 1)) / mean;
    
    double shape;
    if (cv < 0.01) shape = 100.0;
    else shape = 1.2 / cv;
    shape = std::max(0.5, std::min(shape, 20.0));
    
    double log_mean = 0;
    int count = 0;
    for (auto v : shifted) if (v > 0) { log_mean += log(v); count++; }
    log_mean /= count;
    double scale = exp(log_mean + 0.5772 / shape);
    
    if (scale <= 0) scale = mean / 0.885;
    
    return {shape, loc, scale};
}

// ============================================================
// DATA STRUCTURES
// ============================================================

struct DescriptiveStats {
    double mean, median, stddev, variance, skewness, kurtosis;
    double min, max, q1, q3, iqr;
    int n;
    
    void print() const {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  N:          " << n << "\n";
        std::cout << "  Mean:       " << mean << "\n";
        std::cout << "  Median:     " << median << "\n";
        std::cout << "  Std Dev:    " << stddev << "\n";
        std::cout << "  Variance:   " << variance << "\n";
        std::cout << "  Skewness:   " << skewness << "\n";
        std::cout << "  Kurtosis:   " << kurtosis << "\n";
        std::cout << "  Min:        " << min << "\n";
        std::cout << "  Max:        " << max << "\n";
        std::cout << "  Q1:         " << q1 << "\n";
        std::cout << "  Q3:         " << q3 << "\n";
        std::cout << "  IQR:        " << iqr << "\n";
    }
};

struct NormalityResult {
    std::string test_name;
    double statistic;
    double p_value;
    bool is_normal;
};

struct OutlierResult {
    std::vector<int> indices;
    std::vector<double> values;
    std::string method;
    int count;
};

struct ConfidenceInterval {
    std::string parameter;
    double estimate;
    double lower;
    double upper;
    double confidence_level;
};

struct BinomialResult {
    int good_count, bad_count, total_count;
    double defect_proportion;
    double defect_ci_lower, defect_ci_upper;
};

struct CapabilityResult {
    double cpk, cp, ppk, pp, ppm;
    double cpk_ci_lower, cpk_ci_upper;
    bool has_parametric;
};

struct ParameterResult {
    std::string name;
    DescriptiveStats stats;
    std::vector<NormalityResult> normality;
    OutlierResult outliers;
    ConfidenceInterval ci_mean;
    ConfidenceInterval ci_std;
    DistFitResult best_dist;
    BinomialResult binomial;
    CapabilityResult capability;
    std::string method;
};

// ============================================================
// STATISTICAL COMPUTATIONS
// ============================================================

DescriptiveStats compute_descriptive(const std::vector<double>& data) {
    DescriptiveStats s;
    s.n = data.size();
    if (s.n == 0) return s;
    
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    s.mean = std::accumulate(data.begin(), data.end(), 0.0) / s.n;
    s.min = sorted.front();
    s.max = sorted.back();
    s.median = (s.n % 2 == 0) ? (sorted[s.n/2 - 1] + sorted[s.n/2]) / 2.0 : sorted[s.n/2];
    s.q1 = sorted[s.n/4];
    s.q3 = sorted[3*s.n/4];
    s.iqr = s.q3 - s.q1;
    
    double var_sum = 0;
    for (auto v : data) var_sum += (v - s.mean) * (v - s.mean);
    s.variance = var_sum / (s.n - 1);
    s.stddev = sqrt(s.variance);
    
    if (s.n > 2 && s.stddev > 0) {
        double m3 = 0, m4 = 0;
        for (auto v : data) {
            double z = (v - s.mean) / s.stddev;
            m3 += z * z * z;
            m4 += z * z * z * z;
        }
        s.skewness = m3 / s.n;
        s.kurtosis = m4 / s.n - 3.0;
    } else {
        s.skewness = 0;
        s.kurtosis = 0;
    }
    
    return s;
}

std::vector<NormalityResult> normality_tests(const std::vector<double>& data) {
    std::vector<NormalityResult> results;
    
    double sw_p = shapiro_wilk_p(data);
    results.push_back({"Shapiro-Wilk", 0, sw_p, sw_p >= 0.05});
    
    double ad_stat = anderson_darling_stat(data);
    double ad_p = anderson_darling_p(ad_stat);
    results.push_back({"Anderson-Darling", ad_stat, ad_p, ad_p >= 0.05});
    
    return results;
}

OutlierResult detect_outliers_iqr(const std::vector<double>& data) {
    OutlierResult r;
    r.method = "IQR (1.5x)";
    
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    double q1 = sorted[data.size() / 4];
    double q3 = sorted[3 * data.size() / 4];
    double iqr = q3 - q1;
    double lower = q1 - 1.5 * iqr;
    double upper = q3 + 1.5 * iqr;
    
    for (int i = 0; i < (int)data.size(); i++) {
        if (data[i] < lower || data[i] > upper) {
            r.indices.push_back(i);
            r.values.push_back(data[i]);
        }
    }
    r.count = r.indices.size();
    return r;
}

ConfidenceInterval ci_mean_95(const std::vector<double>& data) {
    ConfidenceInterval ci;
    ci.parameter = "Mean";
    ci.confidence_level = 95.0;
    
    DescriptiveStats s = compute_descriptive(data);
    ci.estimate = s.mean;
    
    double se = s.stddev / sqrt(s.n);
    double t_val = t_ppf(0.975, s.n - 1);
    ci.lower = s.mean - t_val * se;
    ci.upper = s.mean + t_val * se;
    
    return ci;
}

ConfidenceInterval ci_std_95(const std::vector<double>& data) {
    ConfidenceInterval ci;
    ci.parameter = "Std Dev";
    ci.confidence_level = 95.0;
    
    DescriptiveStats s = compute_descriptive(data);
    ci.estimate = s.stddev;
    
    double chi2_lower = chi2_ppf(0.025, s.n - 1);
    double chi2_upper = chi2_ppf(0.975, s.n - 1);
    
    ci.lower = sqrt((s.n - 1) * s.variance / chi2_upper);
    ci.upper = sqrt((s.n - 1) * s.variance / chi2_lower);
    
    return ci;
}

// ============================================================
// DISTRIBUTION FITTING (A2-A3)
// ============================================================

DistFitResult find_best_distribution(const std::vector<double>& data) {
    DistFitResult best;
    best.valid = false;
    best.p_value = 0;
    
    struct DistCandidate {
        std::string name;
        std::function<std::function<double(double)>(const std::vector<double>&)> make_cdf;
    };
    
    std::vector<DistCandidate> candidates = {
        {"Normal", 
            [](const std::vector<double>& d) {
                auto [m, s] = fit_normal(d);
                return make_normal_cdf(m, s);
            }
        },
        {"Lognormal",
            [](const std::vector<double>& d) {
                auto [shape, loc, scale] = fit_lognormal(d);
                return make_lognormal_cdf(shape, loc, scale);
            }
        },
        {"Exponential",
            [](const std::vector<double>& d) {
                auto [loc, scale] = fit_exponential(d);
                return make_expon_cdf(loc, scale);
            }
        },
        {"Gamma",
            [](const std::vector<double>& d) {
                auto [shape, loc, scale] = fit_gamma(d);
                return make_gamma_cdf(shape, loc, scale);
            }
        },
        {"Weibull",
            [](const std::vector<double>& d) {
                auto [shape, loc, scale] = fit_weibull(d);
                return make_weibull_cdf(shape, loc, scale);
            }
        }
    };
    
    for (auto& cand : candidates) {
        try {
            auto cdf = cand.make_cdf(data);
            double ks = ks_test_statistic(data, cdf);
            double p = ks_p_value(ks, data.size());
            
            if (p >= 0.05 && p > best.p_value) {
                best.name = cand.name;
                best.p_value = p;
                best.ks_stat = ks;
                best.valid = true;
            }
        } catch (...) {
            continue;
        }
    }
    
    return best;
}

// ============================================================
// BINOMIAL METHOD (A4)
// ============================================================

BinomialResult compute_binomial(const std::vector<double>& data, 
                                 double lsl, double usl) {
    BinomialResult r;
    r.total_count = data.size();
    r.good_count = 0;
    r.bad_count = 0;
    
    for (auto v : data) {
        if (v >= lsl && v <= usl) r.good_count++;
        else r.bad_count++;
    }
    
    r.defect_proportion = (r.total_count > 0) ? (double)r.bad_count / r.total_count : 0;
    
    // Wilson score interval for 95% CI
    double z = 1.96;
    double n = r.total_count;
    double p_hat = r.defect_proportion;
    double denom = 1 + z * z / n;
    double centre = (p_hat + z * z / (2 * n)) / denom;
    double margin = z * sqrt((p_hat * (1 - p_hat) / n) + (z * z / (4 * n * n))) / denom;
    
    r.defect_ci_lower = std::max(0.0, centre - margin);
    r.defect_ci_upper = std::min(1.0, centre + margin);
    
    return r;
}

// ============================================================
// CAPABILITY ANALYSIS (A5)
// ============================================================

CapabilityResult compute_cpk_parametric(const std::vector<double>& data,
                                         double lsl, double usl,
                                         const DistFitResult& dist) {
    CapabilityResult r;
    r.has_parametric = true;
    
    DescriptiveStats s = compute_descriptive(data);
    double mean = s.mean;
    double stddev = s.stddev;
    
    if (stddev == 0) {
        r.cpk = 1e6;
        r.cp = 1e6;
        r.ppm = 0;
        r.cpk_ci_lower = 1e6;
        r.cpk_ci_upper = 1e6;
        return r;
    }
    
    double z_usl = (usl - mean) / stddev;
    double z_lsl = (mean - lsl) / stddev;
    
    r.cp = (usl - lsl) / (6 * stddev);
    r.cpk = std::min(z_usl, z_lsl) / 3.0;
    r.pp = r.cp;
    r.ppk = r.cpk;
    
    // PPM from fitted distribution
    double defect_prob = 0;
    if (dist.valid) {
        std::function<double(double)> cdf;
        if (dist.name == "Normal") {
            cdf = make_normal_cdf(dist.params[0], dist.params[1]);
        } else if (dist.name == "Lognormal") {
            cdf = make_lognormal_cdf(dist.params[0], dist.params[1], dist.params[2]);
        } else if (dist.name == "Exponential") {
            cdf = make_expon_cdf(dist.params[0], dist.params[1]);
        } else if (dist.name == "Gamma") {
            cdf = make_gamma_cdf(dist.params[0], dist.params[1], dist.params[2]);
        } else if (dist.name == "Weibull") {
            cdf = make_weibull_cdf(dist.params[0], dist.params[1], dist.params[2]);
        } else {
            defect_prob = 1.0 - (norm_cdf((usl - mean) / stddev) - norm_cdf((lsl - mean) / stddev));
            goto calc_ppm;
        }
        defect_prob = 1.0 - (cdf(usl) - cdf(lsl));
    } else {
        defect_prob = 1.0 - (norm_cdf((usl - mean) / stddev) - norm_cdf((lsl - mean) / stddev));
    }
    
    calc_ppm:
    r.ppm = defect_prob * 1e6;
    if (r.ppm < 0) r.ppm = 0;
    
    // Cpk confidence interval
    double cpk_var = (r.cpk * r.cpk) / (2 * s.n - 2) + 1.0 / (9 * s.n);
    double cpk_se = sqrt(cpk_var);
    double z = 1.96;
    r.cpk_ci_lower = r.cpk - z * cpk_se;
    r.cpk_ci_upper = r.cpk + z * cpk_se;
    
    return r;
}

CapabilityResult compute_cpk_binomial(const BinomialResult& binom) {
    CapabilityResult r;
    r.has_parametric = false;
    r.cp = 0;
    r.pp = 0;
    r.ppk = 0;
    
    if (binom.defect_proportion == 0) {
        r.cpk = 1.5;
        r.ppm = 0;
        r.cpk_ci_lower = 0;
        r.cpk_ci_upper = 0;
        return r;
    }
    
    double z = norm_ppf(1.0 - binom.defect_proportion / 2.0);
    r.cpk = z / 3.0;
    r.ppm = binom.defect_proportion * 1e6;
    r.cpk_ci_lower = 0;
    r.cpk_ci_upper = 0;
    
    return r;
}

// ============================================================
// ASCII PLOTS
// ============================================================

void print_histogram(const std::vector<double>& data, int bins = 15) {
    if (data.empty()) return;
    
    double min = *std::min_element(data.begin(), data.end());
    double max = *std::max_element(data.begin(), data.end());
    double range = max - min;
    if (range == 0) range = 1;
    double bin_width = range / bins;
    
    std::vector<int> counts(bins, 0);
    for (auto v : data) {
        int idx = std::min((int)((v - min) / bin_width), bins - 1);
        counts[idx]++;
    }
    
    int max_count = *std::max_element(counts.begin(), counts.end());
    int width = 50;
    
    std::cout << "\n  Histogram:\n";
    for (int i = 0; i < bins; i++) {
        double bin_start = min + i * bin_width;
        double bin_end = bin_start + bin_width;
        int bar_len = (max_count > 0) ? (counts[i] * width / max_count) : 0;
        
        std::cout << "  " << std::setw(7) << std::fixed << std::setprecision(2) << bin_start
                  << "-" << std::setw(7) << bin_end << " |";
        for (int j = 0; j < bar_len; j++) std::cout << "█";
        std::cout << " " << counts[i] << "\n";
    }
    std::cout << "  " << std::string(width + 20, '-') << "\n";
}

void print_box_plot(const std::vector<double>& data) {
    if (data.empty()) return;
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    double min = sorted.front();
    double max = sorted.back();
    double q1 = sorted[data.size() / 4];
    double q3 = sorted[3 * data.size() / 4];
    double med = (data.size() % 2 == 0) ? (sorted[data.size()/2 - 1] + sorted[data.size()/2]) / 2.0 : sorted[data.size()/2];
    double iqr = q3 - q1;
    double lf = q1 - 1.5 * iqr;
    double uf = q3 + 1.5 * iqr;
    
    int width = 40;
    double range = max - min;
    if (range == 0) range = 1;
    auto scale = [&](double v) { return (int)((v - min) * width / range); };
    
    std::cout << "\n  Box Plot:\n  ";
    for (int i = 0; i <= width + 4; i++) std::cout << (i == 0 ? '|' : (i == width + 3 ? '|' : '-'));
    
    int p_lf = scale(std::max(min, lf));
    int p_q1 = scale(q1);
    int p_med = scale(med);
    int p_q3 = scale(q3);
    int p_uf = scale(std::min(max, uf));
    
    std::cout << "\n  ";
    for (int i = 0; i <= width + 4; i++) {
        if (i < p_lf || i > p_uf) {
            std::cout << ((i == p_lf || i == p_uf) ? '|' : ' ');
        } else if (i >= p_q1 && i <= p_q3) {
            if (i == p_med) {
                std::cout << (med > (q1 + q3) / 2 ? "┤" : "├");
            } else {
                std::cout << "█";
            }
        } else {
            std::cout << "-";
        }
    }
    
    std::cout << "\n  Min:" << std::setw(7) << min
              << " Q1:" << std::setw(7) << q1
              << " Med:" << std::setw(7) << med
              << " Q3:" << std::setw(7) << q3
              << " Max:" << std::setw(7) << max << "\n";
}

void print_control_chart(const std::vector<double>& data) {
    if (data.size() < 3) return;
    
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    
    std::vector<double> mr;
    for (int i = 1; i < (int)data.size(); i++) {
        mr.push_back(fabs(data[i] - data[i-1]));
    }
    double mr_mean = std::accumulate(mr.begin(), mr.end(), 0.0) / mr.size();
    
    double d2 = 1.128;
    double sigma = mr_mean / d2;
    double ucl = mean + 3 * sigma;
    double lcl = mean - 3 * sigma;
    
    int width = std::min(60, (int)data.size() * 2);
    int height = 12;
    
    double min_val = std::min({lcl * 0.9, *std::min_element(data.begin(), data.end())});
    double max_val = std::max({ucl * 1.1, *std::max_element(data.begin(), data.end())});
    double range = max_val - min_val;
    if (range == 0) range = 1;
    
    auto scale_y = [&](double v) { return (int)((v - min_val) * (height - 1) / range); };
    
    std::cout << "\n  Control Chart (I-MR):\n";
    std::cout << "  UCL = " << std::fixed << std::setprecision(4) << ucl << "\n";
    std::cout << "  CL  = " << mean << "\n";
    std::cout << "  LCL = " << lcl << "\n\n";
    
    for (int row = height - 1; row >= 0; row--) {
        double val = min_val + row * range / (height - 1);
        std::cout << "  " << std::setw(7) << std::fixed << std::setprecision(2) << val << " |";
        
        int row_ucl = scale_y(ucl);
        int row_cl = scale_y(mean);
        int row_lcl = scale_y(lcl);
        
        for (int i = 0; i < std::min((int)data.size(), width); i++) {
            int col = i * width / data.size();
            int row_val = scale_y(data[i]);
            
            if (row_val == row) {
                if (data[i] > ucl || data[i] < lcl) {
                    std::cout << "●";  // Out of control
                } else {
                    std::cout << "●";
                }
            } else if (row == row_ucl && col < data.size()) {
                std::cout << "-";
            } else if (row == row_lcl && col < data.size()) {
                std::cout << "-";
            } else if (row == row_cl && col < data.size()) {
                std::cout << "-";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "         " << std::string(std::min((int)data.size(), width) + 1, '-') << "\n";
}

// ============================================================
// CSV PARSER
// ============================================================

std::vector<std::vector<std::string>> parse_csv(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return data;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            cell.erase(0, cell.find_first_not_of(" \t\r\n"));
            cell.erase(cell.find_last_not_of(" \t\r\n") + 1);
            row.push_back(cell);
        }
        if (!row.empty()) data.push_back(row);
    }
    return data;
}

// ============================================================
// PARAMETER ANALYSIS
// ============================================================

ParameterResult analyze_parameter(const std::string& name,
                                   const std::vector<double>& values,
                                   double lsl, double usl) {
    ParameterResult result;
    result.name = name;
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  Analyzing Parameter: " << name << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "  Spec Limits: LSL=" << lsl << ", USL=" << usl << "\n";
    std::cout << "  Sample Size: " << values.size() << "\n";
    
    result.stats = compute_descriptive(values);
    std::cout << "\n  Descriptive Statistics:\n";
    result.stats.print();
    
    std::cout << "\n  Outlier Detection (IQR):\n";
    result.outliers = detect_outliers_iqr(values);
    std::cout << "    " << result.outliers.count << " outliers found\n";
    if (result.outliers.count > 0) {
        std::cout << "    First 10 indices: ";
        for (int i = 0; i < std::min(10, result.outliers.count); i++) {
            std::cout << result.outliers.indices[i] << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n  Normality Tests:\n";
    result.normality = normality_tests(values);
    for (auto& nt : result.normality) {
        std::cout << "    " << nt.test_name << ": p=" << nt.p_value 
                  << (nt.is_normal ? " (Normal ✓)" : " (Non-normal)") << "\n";
    }
    
    result.ci_mean = ci_mean_95(values);
    result.ci_std = ci_std_95(values);
    std::cout << "\n  Confidence Intervals (95%):\n";
    std::cout << "    Mean: [" << result.ci_mean.lower << ", " << result.ci_mean.upper << "]\n";
    std::cout << "    Std Dev: [" << result.ci_std.lower << ", " << result.ci_std.upper << "]\n";
    
    std::cout << "\n  Distribution Fitting (A2-A3):\n";
    result.best_dist = find_best_distribution(values);
    
    if (result.best_dist.valid) {
        std::cout << "    ✓ Best fit: " << result.best_dist.name 
                  << " (p=" << result.best_dist.p_value << ")\n";
        std::cout << "    KS statistic: " << result.best_dist.ks_stat << "\n";
        result.method = "parametric";
    } else {
        std::cout << "    ✗ No parametric distribution fits (all p<0.05)\n";
        std::cout << "    Using Binomial method (A4)\n";
        result.binomial = compute_binomial(values, lsl, usl);
        std::cout << "    Good: " << result.binomial.good_count 
                  << ", Bad: " << result.binomial.bad_count 
                  << ", Total: " << result.binomial.total_count << "\n";
        std::cout << "    Defect rate: " << (result.binomial.defect_proportion * 100) << "%\n";
        std::cout << "    Defect 95% CI: [" << (result.binomial.defect_ci_lower * 100) 
                  << "%, " << (result.binomial.defect_ci_upper * 100) << "%]\n";
        result.method = "binomial";
    }
    
    std::cout << "\n  Capability Analysis (A5):\n";
    if (result.method == "parametric") {
        result.capability = compute_cpk_parametric(values, lsl, usl, result.best_dist);
    } else {
        result.capability = compute_cpk_binomial(result.binomial);
    }
    
    std::cout << "    Cpk: " << result.capability.cpk << "\n";
    std::cout << "    Cp:  " << result.capability.cp << "\n";
    std::cout << "    Ppk: " << result.capability.ppk << "\n";
    std::cout << "    Pp:  " << result.capability.pp << "\n";
    std::cout << "    PPM: " << (unsigned long long)result.capability.ppm << "\n";
    
    if (result.capability.has_parametric) {
        std::cout << "    Cpk 95% CI: [" << result.capability.cpk_ci_lower 
                  << ", " << result.capability.cpk_ci_upper << "]\n";
    }
    
    print_histogram(values);
    print_box_plot(values);
    print_control_chart(values);
    
    return result;
}

// ============================================================
// MAIN
// ============================================================

int main(int argc, char* argv[]) {
    std::cout << "\n";
    std::cout << "  ╔══════════════════════════════════════════════╗\n";
    std::cout << "  ║       QUALITY ENGINEERING ANALYZER v1.0     ║\n";
    std::cout << "  ║         Mobile Process Capability Suite     ║\n";
    std::cout << "  ╚══════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    std::string filename;
    std::cout << "  Enter CSV filename (or 'sample' for built-in test): ";
    std::getline(std::cin, filename);
    
    filename.erase(0, filename.find_first_not_of(" \t"));
    filename.erase(filename.find_last_not_of(" \t") + 1);
    
    std::map<std::string, std::vector<double>> parameters;
    std::map<std::string, double> lsl_map, usl_map;
    
    if (filename == "sample" || filename.empty()) {
        std::cout << "\n  Using built-in sample dataset (50 rows, 5 params)\n";
        
        std::vector<double> lengths = {
            100.2, 99.8, 101.0, 99.5, 100.8, 100.0, 99.2, 101.5, 100.3, 99.7,
            100.6, 99.4, 100.1, 99.9, 100.4, 99.6, 100.7, 100.0, 99.3, 101.2,
            99.1, 100.9, 100.5, 98.8, 101.1, 101.3, 100.3, 99.8, 100.1, 101.4,
            99.5, 100.0, 100.6, 99.7, 100.4, 99.2, 101.0, 99.6, 100.2, 99.9,
            100.8, 99.4, 100.5, 101.1, 99.8, 100.0, 100.7, 99.3, 100.1, 99.5
        };
        std::vector<double> diameters = {
            5.01, 4.98, 5.10, 5.02, 4.95, 5.00, 5.08, 4.92, 5.04, 4.97,
            5.06, 4.99, 5.03, 4.96, 5.09, 4.94, 5.07, 5.01, 5.05, 4.93,
            5.12, 4.90, 5.00, 5.11, 4.96, 4.99, 5.04, 4.98, 5.02, 4.91,
            5.07, 5.00, 4.97, 5.03, 5.05, 5.10, 4.94, 5.01, 4.99, 5.08,
            4.96, 5.02, 5.06, 4.93, 5.00, 5.05, 4.98, 5.04, 5.01, 4.95
        };
        std::vector<double> hardness = {
            45, 47, 42, 46, 48, 44, 43, 49, 41, 46,
            45, 47, 44, 48, 42, 46, 43, 47, 45, 44,
            48, 41, 46, 43, 47, 44, 46, 42, 48, 44,
            45, 47, 43, 46, 44, 42, 48, 45, 47, 43,
            46, 44, 48, 41, 45, 47, 44, 46, 43, 48
        };
        std::vector<double> temps = {
            22.1, 22.3, 21.9, 22.0, 22.2, 22.1, 21.8, 22.4, 22.0, 22.3,
            22.1, 21.9, 22.2, 22.0, 22.1, 22.3, 21.8, 22.2, 22.0, 22.1,
            21.7, 22.5, 22.2, 21.9, 22.0, 22.4, 22.3, 22.0, 21.8, 22.2,
            22.1, 21.9, 22.3, 22.0, 22.2, 21.7, 22.4, 22.1, 22.0, 22.3,
            21.9, 22.2, 22.0, 22.5, 22.1, 21.8, 22.3, 22.0, 22.2, 22.1
        };
        std::vector<double> pressures = {
            101.3, 101.1, 101.5, 101.2, 100.9, 101.4, 101.0, 101.6, 101.1, 101.3,
            101.2, 101.4, 101.0, 101.5, 101.1, 101.3, 101.2, 101.4, 101.0, 101.6,
            100.8, 101.7, 101.1, 101.2, 101.4, 101.0, 101.3, 101.1, 101.4, 100.9,
            101.2, 101.6, 101.0, 101.5, 101.1, 101.3, 100.8, 101.4, 101.2, 101.0,
            101.5, 101.1, 101.3, 101.6, 101.2, 101.4, 101.0, 101.5, 101.1, 101.3
        };
        
        parameters["length"] = lengths;
        parameters["diameter"] = diameters;
        parameters["hardness"] = hardness;
        parameters["temperature"] = temps;
        parameters["pressure"] = pressures;
        
        lsl_map["length"] = 99.0;
        usl_map["length"] = 101.5;
        lsl_map["diameter"] = 4.85;
        usl_map["diameter"] = 5.15;
        lsl_map["hardness"] = 40.0;
        usl_map["hardness"] = 50.0;
        lsl_map["temperature"] = 21.5;
        usl_map["temperature"] = 22.8;
        lsl_map["pressure"] = 100.5;
        usl_map["pressure"] = 102.0;
        
    } else {
        // Parse actual CSV
        auto csv_data = parse_csv(filename);
        if (csv_data.size() < 2) {
            std::cerr << "Error: CSV must have at least header + 1 data row.\n";
            return 1;
        }
        
        auto headers = csv_data[0];
        for (auto& h : headers) {
            parameters[h] = {};
        }
        
        for (int i = 1; i < (int)csv_data.size(); i++) {
            for (int j = 0; j < (int)csv_data[i].size() && j < (int)headers.size(); j++) {
                try {
                    double val = std::stod(csv_data[i][j]);
                    parameters[headers[j]].push_back(val);
                } catch (...) {
                    // Skip non-numeric
                }
            }
        }
        
        std::cout << "\n  Enter specification limits for each parameter:\n";
        for (auto& [name, vals] : parameters) {
            if (vals.empty()) continue;
            double pmin = *std::min_element(vals.begin(), vals.end());
            double pmax = *std::max_element(vals.begin(), vals.end());
            std::cout << "  " << name << " (n=" << vals.size() 
                      << ", range: " << pmin << " - " << pmax << ")\n";
            
            std::string input;
            std::cout << "    LSL [" << pmin << "]: ";
            std::getline(std::cin, input);
            lsl_map[name] = input.empty() ? pmin : std::stod(input);
            
            std::cout << "    USL [" << pmax << "]: ";
            std::getline(std::cin, input);
            usl_map[name] = input.empty() ? pmax : std::stod(input);
        }
    }
    
    // Run analysis for each parameter
    std::vector<ParameterResult> results;
    for (auto& pair : parameters) {
        const std::string& name = pair.first;
        std::vector<double>& vals = pair.second;
        
        if (vals.size() < 3) {
            std::cout << "  Skipping " << name << " (insufficient data: " << vals.size() << ")\n";
            continue;
        }
        results.push_back(analyze_parameter(name, vals, lsl_map[name], usl_map[name]));
    }
    
    // ============================================================
    // A6: JOINT METRICS
    // ============================================================
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  A6: JOINT METRICS\n";
    std::cout << std::string(60, '=') << "\n";
    
    double joint_cpk = 1e10;
    double joint_ppm_product = 1.0;
    
    for (auto& r : results) {
        if (r.capability.cpk < joint_cpk) joint_cpk = r.capability.cpk;
        joint_ppm_product *= (1.0 - r.capability.ppm / 1e6);
    }
    
    double joint_ppm = (1.0 - joint_ppm_product) * 1e6;
    
    double overall_sigma = 0;
    if (joint_ppm > 0 && joint_ppm < 1e6) {
        double z = norm_ppf(1.0 - joint_ppm / 2e6);
        overall_sigma = z + 1.5;
    } else {
        overall_sigma = joint_cpk * 3.0 + 1.5;
    }
    
    std::cout << "\n  Joint Cpk:      " << std::setprecision(4) << joint_cpk << "\n";
    std::cout << "  Joint PPM:      " << std::setprecision(2) << std::fixed << joint_ppm << "\n";
    std::cout << "  Sigma Level:    " << std::setprecision(3) << overall_sigma << "\n";
    
    std::string quality;
    if (overall_sigma >= 6.0) quality = "World Class (Six Sigma)";
    else if (overall_sigma >= 5.0) quality = "Excellent";
    else if (overall_sigma >= 4.0) quality = "Good (Industry Average)";
    else if (overall_sigma >= 3.0) quality = "Marginal";
    else if (overall_sigma >= 2.0) quality = "Poor";
    else quality = "Critical - Requires Immediate Action";
    
    std::cout << "  Quality Level:  " << quality << "\n";
    
    // ============================================================
    // SUMMARY TABLE
    // ============================================================
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  PARAMETER SUMMARY TABLE\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << std::left;
    std::cout << "  " << std::setw(15) << "Parameter"
              << std::setw(12) << "Method"
              << std::setw(15) << "Distribution"
              << std::setw(10) << "Cpk"
              << std::setw(12) << "PPM"
              << "Status\n";
    std::cout << "  " << std::string(74, '-') << "\n";
    
    for (auto& r : results) {
        std::string status;
        if (r.capability.cpk >= 1.33) status = "Good ✓";
        else if (r.capability.cpk >= 1.0) status = "Marginal ⚠";
        else status = "Poor ✗";
        
        std::cout << "  " << std::setw(15) << r.name
                  << std::setw(12) << r.method
                  << std::setw(15) << (r.best_dist.valid ? r.best_dist.name : "Binomial")
                  << std::setw(10) << std::setprecision(4) << r.capability.cpk
                  << std::setw(12) << std::setprecision(2) << std::fixed << r.capability.ppm
                  << status << "\n";
    }
    
    // ============================================================
    // EXPORT
    // ============================================================
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  REPORT EXPORT\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::string report_name = "qe_report_" + std::to_string(time(nullptr)) + ".txt";
    std::ofstream report(report_name);
    
    if (report.is_open()) {
        time_t now = time(nullptr);
        report << "QUALITY ENGINEERING ANALYSIS REPORT\n";
        report << "Generated: " << ctime(&now);
        report << std::string(60, '=') << "\n\n";
        
        for (auto& r : results) {
            report << "Parameter: " << r.name << "\n";
            report << "  Method: " << r.method << "\n";
            report << "  Distribution: " << (r.best_dist.valid ? r.best_dist.name : "Binomial") << "\n";
            report << "  N: " << r.stats.n << "  Mean: " << r.stats.mean 
                   << "  Std: " << r.stats.stddev << "\n";
            report << "  Cpk: " << r.capability.cpk << "  PPM: " 
                   << (unsigned long long)r.capability.ppm << "\n";
            report << "  LSL: " << lsl_map[r.name] << "  USL: " << usl_map[r.name] << "\n";
            report << "\n";
        }
        
        report << "JOINT METRICS:\n";
        report << "  Joint Cpk: " << joint_cpk << "\n";
        report << "  Joint PPM: " << joint_ppm << "\n";
        report << "  Sigma Level: " << overall_sigma << "\n";
        report << "  Quality: " << quality << "\n";
        
        report.close();
        std::cout << "  ✓ Report saved to: " << report_name << "\n";
        std::cout << "  ✓ Full path: " << report_name << "\n";
    } else {
        std::cout << "  ✗ Could not save report\n";
    }
    
    std::cout << "\n  Press Enter to exit...";
    std::cin.get();
    return 0;
}

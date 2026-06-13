Below is the complete README.md for your GitHub repository. It includes instructions for using both the built‑in sample mode and a separate sample_dataset.csv file that users can create. The entire content is ready to copy and paste.

```markdown
# Quality Engineering Analyzer (QE Analyzer)

**Process Capability & Statistical Analysis Suite for Mobile (Termux)**

This C++ program performs comprehensive quality engineering analysis on CSV data directly in Termux (Android terminal) or any Linux environment. It computes descriptive statistics, normality tests, distribution fitting, process capability (Cpk, PPM), binomial defect analysis, joint metrics, and visual ASCII plots – all without external dependencies.

## Features

- 📊 **Descriptive Statistics** – mean, median, variance, skewness, kurtosis, quartiles, IQR
- 🔬 **Normality Tests** – Shapiro‑Wilk, Anderson‑Darling
- 📈 **Distribution Fitting** – Normal, Lognormal, Exponential, Gamma, Weibull (KS test based)
- ✅ **Capability Analysis** – Cpk, Cp, Ppk, Pp, PPM with 95% confidence intervals
- 🎯 **Binomial Method** – for non‑normal data (defect proportion and Wilson score CI)
- 🔍 **Outlier Detection** – IQR (1.5×) method with indices listing
- 📉 **ASCII Plots** – Histogram, Box Plot, I‑MR Control Chart (in terminal)
- 📑 **Joint Metrics** – Overall Cpk, Sigma Level, Quality classification
- 📁 **CSV Input** – flexible header recognition, skip non‑numeric cells
- 📄 **Auto‑generated Report** – plain text export with timestamp

## Prerequisites (Termux)

Install the required packages:

```bash
pkg update && pkg upgrade
pkg install git build-essential
```

These provide g++ (C++17 compiler) and standard build tools.

Installation & Compilation

Clone the repository and compile the source:

```bash
git clone https://github.com/brhanumehari/quality_analyzer.git
cd quality_analyzer
g++ -std=c++17 -O3 -o quality_analyzer quality_analyzer.cpp -lm
```

Note: -lm links the math library. The compiler must support C++17.

Usage

Run the program:

```bash
./quality_analyzer
```

You will be prompted to enter a CSV file path. You can type sample to run a built‑in demo dataset (50 rows, 5 parameters), or provide the path to your own CSV file.

Using the Built‑in Sample (no file needed)

Just type sample at the prompt. The program will generate internal data and analyse five parameters (length, diameter, hardness, temperature, pressure) with predefined specification limits.

Using Your Own CSV File

If you provide a CSV file, the program will ask for specification limits (LSL, USL) for each numeric column. Defaults are taken from the data range.

Input CSV Format

· First row must contain column headers.
· Following rows contain numeric values (non‑numeric cells are ignored).
· Columns are separated by commas (,).

Example: sample_dataset.csv

You can create a sample CSV file yourself. Copy the content below into a file named sample_dataset.csv:

```csv
length,diameter,hardness,temperature,pressure
100.2,5.01,45,22.1,101.3
99.8,4.98,47,22.3,101.1
101.0,5.10,42,21.9,101.5
99.5,5.02,46,22.0,101.2
100.8,4.95,48,22.2,100.9
100.0,5.00,44,22.1,101.4
99.2,5.08,43,21.8,101.0
101.5,4.92,49,22.4,101.6
100.3,5.04,41,22.0,101.1
99.7,4.97,46,22.3,101.3
100.6,5.06,45,22.1,101.2
99.4,4.99,47,21.9,101.4
100.1,5.03,44,22.2,101.0
99.9,4.96,48,22.0,101.5
100.4,5.09,42,22.1,101.1
99.6,4.94,46,22.3,100.9
100.7,5.07,43,21.8,101.2
100.0,5.01,47,22.2,101.4
99.3,5.05,45,22.0,101.0
101.2,4.93,44,22.1,101.6
99.1,5.12,48,21.7,100.8
100.9,4.90,41,22.5,101.7
100.5,5.00,46,22.2,101.1
98.8,5.11,43,21.9,101.2
101.1,4.96,47,22.0,101.4
101.3,4.99,44,22.4,101.0
100.3,5.04,46,22.3,101.3
99.8,4.98,42,22.0,101.1
100.1,5.02,48,21.8,101.4
101.4,4.91,41,22.2,100.9
99.5,5.07,45,22.1,101.2
100.0,5.00,47,21.9,101.6
100.6,4.97,43,22.3,101.0
99.7,5.03,46,22.0,101.5
100.4,5.05,44,22.2,101.1
99.2,5.10,42,21.7,101.3
101.0,4.94,48,22.4,100.8
99.6,5.01,45,22.1,101.4
100.2,4.99,47,22.0,101.2
99.9,5.08,43,22.3,101.0
100.8,4.96,46,21.9,101.5
99.4,5.02,44,22.2,101.1
100.5,5.06,48,22.0,101.3
101.1,4.93,41,22.5,101.6
99.8,5.00,45,22.1,101.2
100.0,5.05,47,21.8,101.4
100.7,4.98,44,22.3,101.0
99.3,5.04,46,22.0,101.5
100.1,5.01,43,22.2,101.1
99.5,4.95,48,22.1,101.3
```

Place this file in a convenient location (e.g., ~/storage/downloads/sample_dataset.csv on Termux after running termux-setup-storage). Then run:

```bash
./quality_analyzer
# At the prompt, enter: /data/data/com.termux/files/home/storage/downloads/sample_dataset.csv
```

The program will ask you for LSL and USL for each parameter. Suggested limits:

Parameter LSL USL
length 99.0 101.5
diameter 4.85 5.15
hardness 40 50
temperature 21.5 22.8
pressure 100.5 102.0

Sample Run (with built‑in sample command)

```bash
./quality_analyzer
Enter CSV filename (or 'sample'): sample
```

The program will analyse five parameters and produce:

· Descriptive stats
· Normality test results
· Best fitting distribution (if any)
· Capability indices (Cpk, PPM)
· ASCII histogram, box plot, control chart
· Final joint metrics and a summary table
· A report file (qe_report_<timestamp>.txt)

Output Interpretation

Capability Indices

Cpk Value Interpretation
≥ 1.33 Good (capable)
1.00 – 1.32 Marginal
< 1.00 Poor (incapable)

Sigma Level (from joint metrics)

Sigma Quality Level
≥ 6.0 World Class (Six Sigma)
5.0 – 5.9 Excellent
4.0 – 4.9 Good (Industry Average)
3.0 – 3.9 Marginal
2.0 – 2.9 Poor
< 2.0 Critical – Immediate Action

ASCII Control Chart

· ● – individual measurement
· - – control limits (UCL / CL / LCL)
· ● marks an out‑of‑control point

Report Generation

After analysis, a text report is automatically saved as qe_report_<timestamp>.txt in the current directory. It contains parameter‑wise summary, joint Cpk, sigma level, and quality classification.

Termux Quick Start

1. Install Termux from F‑Droid (not the Play Store).
2. Grant storage access:
   ```bash
   termux-setup-storage
   ```
3. Install compiler and git:
   ```bash
   pkg update && pkg install git build-essential
   ```
4. Clone, compile, and run as described above.
5. For easy access, copy your CSV files into ~/storage/downloads/.

Uninstalling

Simply delete the compiled binary and the source folder:

```bash
rm -f quality_analyzer
cd .. && rm -rf quality_analyzer
```

Troubleshooting

Compilation Errors

· Make sure you are using Termux and have installed build-essential.
· If you get undefined reference to 'std::...', try:
  ```bash
  pkg install libc++ clang
  clang++ -std=c++17 -O3 -o quality_analyzer quality_analyzer.cpp -lm
  ```

CSV Not Found

Provide the full path, e.g., /data/data/com.termux/files/home/storage/downloads/mydata.csv or use ~/storage/downloads/mydata.csv.

Non‑numeric Data

The program skips cells that cannot be converted to double. If a column contains only text, it will be ignored.

Contributing

Feel free to open issues or pull requests on GitHub. Possible improvements:

· Parallel processing for large datasets
· More distribution families (Beta, Logistic)
· HTML/PDF report export
· Direct reading from Google Sheets or SQLite

License

MIT License – use freely for academic, industrial, or personal quality analysis.

Author

Brhanu Mehari
GitHub: @brhanumehari
 

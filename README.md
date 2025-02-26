# VPD-thresholds_paper_codes

Codes for paper *"Global vegetation responses to elevated CO2 vary across VPD thresholds"*

## 1. System Requirements
- **Software Dependencies & Operating Systems:**
  - R
  - RStudio (optional, recommended)
  - Required main R packages:
    - `ggplot2` (for visualization)
    - `grf` (for causal analysis)
    - `h2o` (for machine learning)
- **Versions the codes Has Been Tested On:**
  - R 4.4.2
- **Required Non-Standard Hardware:**
  - No specialized hardware required. Code runs on any standard desktop/laptop.

## 2. Demo
- **Instructions to Run on Data:**
  1. Open `R` or `RStudio`.
  2. Set the working directory where your code and data are located:
     ```r
     setwd("/path/to/your/project")
     ```
  3. Load and preprocess the dataset:
     ```r
     load("demo_data_monthly.RData")  # Runs data preprocessing
     ```
  4. Run the main analysis script:
    
- **Expected Output:**
  - Estimated eCO2 effects on NIRv
  - ML models outputs
  - VPD thresholds
  - Estimated Average Treatment Effects (ATE)
  - VPD sensitivity
  - VPD sensitivity trend analysis
  - Visualizations

## 3. Instructions for Use
- **How to Run the Software on Your Data:**
  - Replace the dataset in `data/` directory with your own dataset (ensure it follows the expected format).
  - Update paths in the scripts if needed.
  - Run the scripts in sequence:

- **(OPTIONAL) Reproduction Instructions:**
  - If you wish to reproduce the exact results in the manuscript, use the provided dataset `demo_data` and execute the steps above without modifications.
  - For ensuring reproducibility, we recommend setting a random seed:
    ```r
    set.seed(1234)
    ```
  - If using parallel computing, ensure you have `doParallel` installed and specify the number of cores:
    ```r
    library(doParallel)
    cl <- makeCluster(detectCores() - 1)  
    registerDoParallel(cl)
    ```

---

**We encourage you to include instructions for reproducing all the quantitative results in the manuscript.**

# VPD-thresholds_paper_codes

Codes for paper *"Global vegetation responses to elevated CO2 vary across VPD thresholds"*

---

## 1. System Requirements

* **Software Dependencies & Operating Systems:**

  * R (≥ 4.4.2)
  * RStudio (recommended)
  * Required R packages:

    * `ggplot2` (visualization)
    * `grf` (causal analysis)
    * `h2o` (machine learning)
    * `doParallel` (optional, for parallel computing)

* **Tested Environment:**

  * R 4.4.2 on Windows/macOS/Linux

* **Hardware Requirements:**

  * No specialized hardware required. Runs on standard desktop/laptop.

---

## 2. Project Structure

This repository is organized as follows:

* `*.Rproj` — RStudio project file (recommended entry point)
* `*.R` — analysis scripts

Using the `.Rproj` file ensures all paths are handled **relative to the project root**, avoiding the need for absolute paths.

---

## 3. Demo

### How to Run

1. **Open the project**

   * Double-click the `.Rproj` file
   * Or open it via RStudio (`File → Open Project`)

2. **Load and preprocess data**

   ```r
   load("suppl_data/data_monthly.RData")
   ```

3. **Run the main analysis scripts**
   (run scripts in order, e.g.)

---

### Expected Output

* Estimated eCO2 effects
* Machine learning model outputs
* VPD thresholds
* Estimated Average Treatment Effects (ATE)
* VPD sensitivity
* VPD sensitivity trend analysis
* Visualizations

---

### Expected Run Time

* Data loading: ~30 seconds
* Model training: depends on CPU and RAM
* Visualization: ~10 seconds


---

## Notes

* All paths are relative to the project root

---


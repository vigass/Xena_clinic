## Code Explanation

### 1. Setup
- **Clear workspace and set working directory**: 
  The code first removes all objects from the current workspace using `rm(list = ls())` and sets the working directory to `"D:/Project/TCGA/ESCA/clinic"` using `setwd()`.
  
- **Load necessary libraries**: 
  The required libraries `data.table` and `tidyverse` are loaded. The `data.table` library is used for fast data manipulation, while `tidyverse` provides a set of tools for data analysis and manipulation.

### 2. Load Data
- **Clinical data**: 
  The clinical data is loaded from a `.tsv.gz` file (`TCGA-ESCA.GDC_phenotype.tsv.gz`) using `fread()`.
  
- **Survival data**: 
  The survival data is loaded from a `.tsv` file (`TCGA-ESCA.survival.tsv`) using `read.delim()`.

### 3. Preliminary Exploration
- **Column names**: 
  The column names of the clinical data are stored in the `tmp` data frame.
  
- **Stage, M, T, N, and Grade Distribution**: 
  Tables are generated using `table()` to examine the distributions of the tumor stage, pathologic M, T, and N stages, as well as the tumor grade.

### 4. Data Selection and Column Renaming
- **Select relevant columns**: 
  A subset of the clinical data is created by selecting the following columns:
  - `submitter_id.samples`: Sample ID
  - `neoplasm_histologic_grade`: Tumor grade
  - `pathologic_T`: T stage
  - `pathologic_M`: M stage
  - `pathologic_N`: N stage
  - `age_at_index.demographic`: Age
  - `gender.demographic`: Gender
  - `tumor_stage.diagnoses`: Tumor stage
  
- **Rename columns**: 
  The selected columns are renamed for clarity:
  - `ID`, `grade`, `T_stage`, `M_stage`, `N_stage`, `age`, `gender`, `stage`

### 5. Data Cleaning
- **Grade cleaning**: 
  The `grade` column is filtered to remove rows where `grade` is `"GX"`.

- **Stage cleaning**: 
  The `stage` column is filtered to remove values such as `"not reported"` and `"stage x"`. Stages are cleaned by removing extra characters like `"abc"` and converting Roman numerals to numeric values based on a mapping:
  - Stage I -> 1
  - Stage II -> 2
  - Stage III -> 3
  - Stage IV -> 4

- **T stage cleaning**: 
  The `T_stage` column is cleaned by removing any letters (like "abcd") after the numeric stage.

- **M stage cleaning**: 
  The `M_stage` column is cleaned by removing any text within parentheses and changing `"M1a"` to `"M1"`. Only rows where the `M_stage` is either `"M0"` or `"M1"` are retained.

- **N stage cleaning**: 
  The `N_stage` column is cleaned by removing any text within parentheses and any extra letters (like "abcmi"). Rows with `N_stage` equal to `"NX"` are filtered out.

### 6. Merge Clinical and Survival Data
- **Filter survival data**: 
  The survival data is filtered to only include samples that are present in the cleaned `meta` data.
  
- **Merge data**: 
  The filtered survival data is merged with the cleaned clinical data (`meta`) using the `sample` column from the survival data and the `ID` column from the clinical data.

### 7. Save Processed Data
- The final processed data (`pd`) is saved as an R data file (`pd.rdata`) using `save()`.


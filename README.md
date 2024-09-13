## Code Explanation

### 1. Setup
- **Clear workspace and set working directory**: 
  The code first removes all objects from the current workspace using `rm(list = ls())` and sets the working directory to `"D:/Project/TCGA/ESCA/clinic"` using `setwd()`.
  
- **Load necessary libraries**: 
  The required libraries `data.table` and `tidyverse` are loaded. The `data.table` library is used for fast data manipulation, while `tidyverse` provides a set of tools for data analysis and manipulation.

### 2. Load Data
- **Clinical data**: 
  The clinical data is loaded from a `.tsv.gz` file (`TCGA-ESCA.GDC_phenotype.tsv.gz`) using `fread()`.

  `TCGA-ESCA.GDC_phenotype.tsv.gz`:download by (https://xenabrowser.net/datapages/?dataset=TCGA-ESCA.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
  
- **Survival data**: 
  The survival data is loaded from a `.tsv` file (`TCGA-ESCA.survival.tsv`) using `read.delim()`.

  `TCGA-ESCA.survival.tsv`: download by (https://xenabrowser.net/datapages/?dataset=TCGA-ESCA.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

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


## 代码解释

### 1. 设置
- **清空工作空间并设置工作目录**：
  使用 `rm(list = ls())` 清除当前工作空间中的所有对象，并使用 `setwd()` 将工作目录设置为 `"D:/Project/TCGA/ESCA/clinic"`。
  
- **加载必要的库**：
  加载了 `data.table` 和 `tidyverse` 库。`data.table` 用于快速的数据处理，而 `tidyverse` 提供了一套用于数据分析和处理的工具。

### 2. 加载数据
- **临床数据**：
  使用 `fread()` 从 `.tsv.gz` 文件 (`TCGA-ESCA.GDC_phenotype.tsv.gz`) 中加载临床数据。
  
- **生存数据**：
  使用 `read.delim()` 从 `.tsv` 文件 (`TCGA-ESCA.survival.tsv`) 中加载生存数据。

### 3. 初步探索
- **列名**：
  临床数据的列名存储在 `tmp` 数据框中。
  
- **分期、M、T、N 和肿瘤分级分布**：
  使用 `table()` 函数生成表格，查看肿瘤分期、病理M、T、N分期以及肿瘤分级的分布情况。

### 4. 数据选择和列重命名
- **选择相关列**：
  从临床数据中选择以下列生成子集：
  - `submitter_id.samples`: 样本ID
  - `neoplasm_histologic_grade`: 肿瘤分级
  - `pathologic_T`: T分期
  - `pathologic_M`: M分期
  - `pathologic_N`: N分期
  - `age_at_index.demographic`: 年龄
  - `gender.demographic`: 性别
  - `tumor_stage.diagnoses`: 肿瘤分期
  
- **重命名列**：
  为清晰起见，将选定的列重命名为：
  - `ID`, `grade`, `T_stage`, `M_stage`, `N_stage`, `age`, `gender`, `stage`

### 5. 数据清理
- **清理分级**：
  使用 `filter()` 函数过滤掉 `grade` 列中为 `"GX"` 的行。

- **清理分期**：
  过滤掉 `stage` 列中值为 `"not reported"` 和 `"stage x"` 的行。然后通过删除 `"abc"` 等字符并根据映射将罗马数字转换为阿拉伯数字：
  - I期 -> 1
  - II期 -> 2
  - III期 -> 3
  - IV期 -> 4

- **清理T分期**：
  使用 `gsub()` 删除 `T_stage` 列中的任何字母（如 "abcd"）后面的字符。

- **清理M分期**：
  使用 `gsub()` 删除 `M_stage` 列中的括号内容，并将 `"M1a"` 改为 `"M1"`。仅保留 `M_stage` 为 `"M0"` 或 `"M1"` 的行。

- **清理N分期**：
  使用 `gsub()` 删除 `N_stage` 列中的括号内容和任何多余字母（如 "abcmi"）。过滤掉 `N_stage` 等于 `"NX"` 的行。

### 6. 合并临床和生存数据
- **筛选生存数据**：
  筛选生存数据，保留样本ID在清理后的 `meta` 数据中的样本。
  
- **合并数据**：
  使用 `merge()` 函数，通过生存数据的 `sample` 列和临床数据的 `ID` 列，合并筛选后的生存数据和清理后的临床数据。

### 7. 保存处理后的数据
- 使用 `save()` 将最终处理好的数据（`pd`）保存为 R 数据文件 (`pd.rdata`)。



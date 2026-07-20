# Data Dictionary

## 1. Purpose

This document describes the datasets used in the reproducibility materials for the paper "A general randomized test for alpha" written by
Daniele Massacci, Daniele Massacci, Lucio Sarno, Lorenzo Trapani, and Pierluigi Vallarino.

## 2. Data files

### Dataset: `Asset returns`

**Files:**
- `Empirics/data_S&P500_fullSample_min5Yrs.RData`

**Source and access:**  
Retrieved from Datastream using the license of Erasmus University Rotterdam

**Description:**  
Monthly gross returns (in percentage points, e.g. a return of 5% is stored as 5)
on constituents of the S&P500  between January 1980 and December 2024. 
Returns are continuous data.

Loading the dataset returns a list of two objects called "lData":

  - The first object is called "dates" and is a vector of R Dates formatted as "yyyy-mm-dd";
  - The second object is called "ret" and is a (539, 1019) dimensional matrix containing the gross returns of interest. Column names are the names of the
    firms while rows represent dates. Missing returns are labelled by "NA".

### Dataset: `Fama French factors and risk free rates`

**Files:**
- `Empirics/FF_Data_new.csv`

**Source and access:**  
Retrieved from the website of Professor Kenneth French  (https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)

**Description:**  
Monthly values of the Fama French factors and the risk-free rate used in our analysis.
Loading this dataset returns a list of eight elements. 

  - The first object is called "X" and is a vector of integers representing dates as yyyymmdd;
  - The second one is called "Mkt.RF" and contains market excess returns in percentage points (continuous variable);
  - The third one is called "SMB" and contains values of the size factor in percentage points (continuous variable);
  - The fourth one is called "HML" and contains values of the value factor in percentage points (continuous variable);
  - The fifth one is called "RMW" and contains values of the profitability factor in percentage points (continuous variable);
  - The sixth one is called "CMA" and contains values of the investment strategy factor in percentage points (continuous variable);
  - The seventh one is called "RF" and contains values of the risk-free rate in percentage points (continuous variable);
  - The eighth one is called "MOM" and contains values of the momentum factor in percentage points (continuous variable);
    
### Variable transformations

The empirical analysis uses the Fama French factors as they are (and hence, as downloaded from the relevant sources). Since the factor models that we test
use excess returns as dependent variables, the risk-free rate is subtracted from  the gross returns imported from the .Rdata file.



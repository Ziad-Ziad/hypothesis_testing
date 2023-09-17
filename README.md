# Hypothesis Testing

The module supports a wide range of 1 and 2-sample statistical tests (A/B tests), allowing you to analyze and draw insights from data effectively. It includes the following test statistics:
- Z-test
- Independent T-test
- Proportion Test
- Paired T-test
- F-test
- Chi-squared Contingency Test

## Usage

To use this script, follow these steps:

1. Ensure that you have the necessary input data files:

   - For Z-test, Independent T-test, Paired T-test, and F-test, you'll need two data files: `control.txt` and `variation.txt`.
   - For the Chi-squared Contingency Test, you'll need a data file `table.txt`.
   - For the Proportion Test, you'll need to specify the data as separate lines in a file named `data.txt`, where each line contains the values for obs1, n1, obs2, and n2 in that order.

2. Run the script using the following command:
```shell
python main.py <test_type> [--alpha=<alpha>] [--alter=<alter>]
```
- `<test_type>`: Choose one of the following options:
  - `z_test`: Z-test for two samples.
  - `t_independent`: Independent two-sample T-test.
  - `proportion_test`: Two-sample proportion test.
  - `paired_test`: Paired T-test.
  - `f_test`: F-test for comparing variances.
  - `chi2_conting`: Chi-squared contingency test.
- `--alpha <alpha>`: Set the significance level (default is 0.05).
- `--alter <alter>`: Specify the alternative hypothesis (options: `two-tailed`, `right`, `left`). Note that not all tests support all alternative hypotheses.

3. View the test results, which will include the test statistic, p-value, chosen alpha level, and alternative hypothesis.

## Input Data Files

- `control.txt`: Data file for the control group.
- `variation.txt`: Data file for the variation group.
- `table.txt`: Data file for the Chi-squared contingency test.
- `data.txt`: Data file for the Proportion Test (obs1, n1, obs2, n2 on separate lines).

## Requirements
To run this script, you'll need the following dependencies:

- Python 3.x
- NumPy
- SciPy

## Example

Here's an example of how to run the script for a Z-test:

```shell
python main.py z_test --alpha=0.01 --alter='right'
```

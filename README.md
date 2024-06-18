# RNN-THE
# Datasets
In fold "Testing datasets" where "Testing data examples" can be directly used for experimental verifying (Physionet Challenge 2019, by default should be located at "D:\") and "Testing datasets" contains all four testing longitudinal patient datasets (need to be randomly split into multiple training and testing sets before experimental verifying)
# Usage
Directly run "testingCommand.bat"
"./pseudo-HE and nonsecure SRNN.exe" 9 0
"./pseudo-HE and nonsecure RNN.exe" 9 0)
The first number indicates the number of training-testing sets (by default 9 means 10 training-testing sets), the second number indicates the hyperparameter of l2-regularization. Then, four files named like "Convergence_RNN0_0.000000.csv", "Convergence_SRNN0_0.000000.csv", "Result_RNN0_0.000000.csv", "Result_SRNN0_0.000000.csv" will appear at "D:\".
The former two indicates the variance of the training and testing sets' overall loss with epoch time grows, the latter two indicates the predictive value of testing sets with epoch time grows, which can be used to evaluate prediction performance.
To change the version between pseudo-HE and nonsecure one, just change the activation function and regenerate the binary files.
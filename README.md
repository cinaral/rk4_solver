# ```rk4_solver```: Runge-Kutta 4th Order Solver
Runge-Kutta 4th Order Method ODE Solver with events. This is a header-only library.

## Testing
Reference data is needed for the tests. By default they can be found in ```test/dat/```. 

You may need to generate new reference data in order to update the existing tests or to add new tests. ```test/matlab/run_all_tests.m``` will generate reference data in ```dat/``` if you have access to MATLAB. Then the generated data (```*.dat```) can be copied into ```test/dat/```. 

The ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```.
```MATLAB
writematrix(matrix, file, 'Delimiter', ',');  
matrix = readmatrix(file);  
```
 Each row in file corresponds to a matrix row. For example, if the test reference data is an ```N``` by ```M``` matrix then the ```*.dat``` file should contain ```N``` rows and ```M``` numbers at each row.
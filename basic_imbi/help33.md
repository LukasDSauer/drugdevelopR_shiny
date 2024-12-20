#### Values presented in table
* expected utility (u)
* assumed true standardized difference in means (&Delta;)
* threshold value for the decision rule to go to phase III (&kappa;)
* total sample size for phase II (n<sub>2</sub>)
* total expected sample size for phase III (n<sub>3</sub>, rounded to the next equal natural number)
* total expected sample size in the program (n = n<sub>2</sub> + n<sub>3</sub>)
* maximal costs of the program (K)
* probability to go to phase III (p<sub>go</sub>)
* probability of a successful program (sProg)
* probability of a successful program with small treatment effect in Phase III (sProg1)
* probability of a successful program with medium treatment effect in Phase III (sProg2)
* probability of a successful program with large treatment effect in Phase III (sProg3)
* expected costs for phase II (K2)
* expected costs for phase III (K3)

and further input parameters.


#### References

Cohen, J. (1988). Statistical power analysis for the behavioral sciences.

IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at https://www.iqwig.de/de/methoden/methodenpapier.3020.html, assessed last 15.05.19.


#### Additional Features
The true underlying standardized difference in means can also be modelled by assuming a prior distribution (https://web.imbi.uni-heidelberg.de/prior/). Optimal planning can then be done with the function optimal_normal() of the R package drugdevelopR available at: https://github.com/Sterniii3/drugdevelopR.

#### Note

If the server is busy, you may need to double click the "Go"-button in order to see the updated plot.

#### Maintainer

Stella Erdmann, Institute of Medical Biometry, University of Heidelberg, email: erdmann@imbi.uni-heidelberg.de.

Version 1.0




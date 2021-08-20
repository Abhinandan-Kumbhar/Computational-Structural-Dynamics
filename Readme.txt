There are total 4 octave codes.

For solving 3 DOF example for validation.
1. GMAM for example problem

For solving project problem:
1.Newmarkforproject - Octave code to find required response at 1891 node using Newmark beta method.
2. GMAMprojectproblem - Octave code to find required response at 1891 node using GMAM.
3. Overall   - Combined Octave code for Newmark method and GMAM

Run GMAMprojectproblem code for results of GMAM applied to project problem.
Run Newmarkforproject code for results Newmark method applied to project problem.

Instead of running above two codes separately, 
Run 'Overall' code, it is used to find comprehensive results shown in report.
Though time consuming (15-20 min) it helps to get insight about the accuracy of GMAM.

[   1,   1]: 7.528e+09 [   1,   2]:-2.817e+08 [   1,   3]: 3.722e+04
[   2,   2]: 7.006e+09 [   2,   3]:-6.918e+04 [   2,   4]:-3.168e+08 
 
 

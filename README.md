## PFG ##

PFG is a parallel formula generator to calculate the possible elemental compositions for a given mass. PFG is implemented in C++. It can be compiled by MSVC and GCC easily and run smoothly in both Windows and Linux.

## Usage ##
```
-h or --help       The help screen.
-m mass            Set mass.
-t tol             Set tolerance to tol 'ppm' ( default 5 ).
-u unit            Set unit of tol, could be 'Da' or 'ppm' (default ppm).
-c charge          Set charge to be calculated.
-r rules           Set rules to constrain formulas
-f file            Set filename to stored the generated formulas (default is stdout).
--X a-b            Set atom range a (min) to b (max) of element X.
                   some of the valid elements:
          X           key      mass(6 decimals shown)
     ----------------------------------------------------
          C           --C           12.000000
          13C         --13C         13.003355
          H           --H           1.007825
          D           --D           2.014102
          N           --N           14.003074
          15N         --15N         15.000109
          O           --O           15.994915
          P           --P           30.973762
          S           --S           31.972071
          F           --F           18.998403
          Cl          --Cl          34.968853
          Br          --Br          78.918338
          I           --I           126.904468
          Si          --Si          27.976927
          Na          --Na          22.989770
          K           --K           38.963707
--agentformula af  Set agent formula.
--agentcharge  ac  Set the charge of the agent formula. 
```
## Example ##
```
PFG  'Mesuximide' -m 203.0946 -t 5 --C 0-20 --H 0-40 --O 0-5 --N 0-5 -f test.txt
```
the result were stored in the test.txt

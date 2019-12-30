# wave3

**wave3** is a small program used in Beleza S, Santos AM, McEvoy B, Alves I, Martinho C, Cameron E, Shriver MD, Parra EJ, Rocha J (2013). The timing of pigmentation lightening in Europeans. *Molecular Biology and Evolution*, **30(1):** 24-35.

## Running 

```bash
$ wave3 <configuration_file> <results_file>

```

It produces a resulting text (ASCII) file with nine columns, and a line for each simulation that succeeded. The columns represent

1) *s_number* total number of haplotypes carrying the mutation S
2) *different_s* number of different mutant S haplotypes
3) *original_s* number of haplotypes similar to the original one (11111111S11111111)
4) *heter* heterozigozity
5) *p* p, or the frequency of original mutant haplotyypes within mutant haplotypes
6) *size* population size
7) *sel_coef* selection coeficient
8) *g* generation when the simulation stopped (all stop criteria reached)
9) *s_freq* frequency of the mutant (S) haplotype

```ascii
s_number;different_s;original_s;heter;p;size;sel_coef;g;s_freq
1477;8;73;0.724853;0.0494245;2756;0.0124275;398;0.535922
1058;10;91;0.728926;0.0860113;2284;0.0319153;143;0.463222
1533;8;83;0.653928;0.0541422;2553;0.0327738;284;0.60047
1376;10;163;0.730048;0.118459;2676;0.0148885;352;0.5142
...
```

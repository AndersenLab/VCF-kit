# Overview

`tajima` can be used to calculate [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D). For an explanation of what Tajima's D is, see [this excellent video](https://www.youtube.com/watch?v=wiyay4YMq2A) by Mohamad Noor.

# Usage

Generate a Tajima's D statistic using a sliding window or across bins. 

__Parameters__:

* __window-size__ - Size of window from in which to calculate Tajima's D. 
* __step-size__ - Size of step taken.
* __--sliding__ - Fluidly slide along genome, capturing every window of a given `window-size`. Equivelent to `step-size` = 1;
* __--no-header__ - Outputs results without a header. 
* __--extra__ - Adds on `filename`, `window-size` and `step-size` as additional columns. Useful for comparing different files / parameters. 

![Tajima](tajima.png)

<small>The figure above illustrates the types of windows over which Tajima's D can be calculated. 

### Bin calculation

If you set the `window-size` and `step-size` as the same value, the bins will not overlap.

```
vk tajima 1000 1000 <vcf>
```

The code above will calculate Tajima's D using 100,000 bp bins across the genome.

| CHROM   |   BIN_START |   BIN_END |   N_Sites |   N_SNPs |    TajimaD |
|:--------|------------:|----------:|----------:|---------:|-----------:|
| I       |           0 |    100000 |         1 |        1 | -1.59468   |
| I       |      100000 |    200000 |         2 |        2 | -0.0291034 |
| I       |      200000 |    300000 |         2 |        2 | -1.53214   |
| I       |      300000 |    400000 |         2 |        2 | -1.53214   |
| I       |      400000 |    500000 |         2 |        2 | -0.670611  |
| I       |      500000 |    600000 |         2 |        1 | -1.59468   |



### Sliding window

When run, the code below will calculate Tajima's D across a 100,000 bp sliding window that moves 1,000 bp with each iteratino.

```
vk tajima 100000 1000 <vcf>
```

| CHROM   |   BIN_START |   BIN_END |   N_Sites |   N_SNPs |   TajimaD |
|:--------|------------:|----------:|----------:|---------:|----------:|
| I       |           1 |    100000 |         1 |        1 |  -1.59468 |
| I       |        1001 |    101000 |         1 |        1 |  -1.59468 |
| I       |        2001 |    102000 |         1 |        1 |  -1.59468 |
| I       |        3001 |    103000 |         1 |        1 |  -1.59468 |
| I       |        4001 |    104000 |         1 |        1 |  -1.59468 |
| I       |        5001 |    105000 |         1 |        1 |  -1.59468 |


### Continous sliding window

When run, the code below will calculate Tajima's D across a 100,000 bp sliding window that captures every unique bin of variants that fall within 100,000 bp of one another.

```
vk tajima 100000 --sliding <vcf>
```

| CHROM   |   BIN_START |   BIN_END |   N_Sites |   N_SNPs |    TajimaD |
|:--------|------------:|----------:|----------:|---------:|-----------:|
| I       |           0 |    100000 |         1 |        1 | -1.59468   |
| I       |           0 |    100000 |         2 |        2 | -1.53214   |
| I       |       90777 |    190777 |         2 |        2 | -0.0291034 |
| I       |      154576 |    254576 |         2 |        2 | -0.589097  |
| I       |      207871 |    307871 |         2 |        2 | -1.53214   |



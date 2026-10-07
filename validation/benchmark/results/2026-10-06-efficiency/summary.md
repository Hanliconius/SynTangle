# Efficiency overview

| Budget | Start | Variant | Complete/expected | Proven | Median wall s | Median LNS s | Median global s |
|---:|---|---|---:|---:|---:|---:|---:|
| 150 | pre_hybrid | equal_nohint | 18/18 | 12 | 2.96 | 0.00 | 0.75 |
| 150 | pre_hybrid | weighted_adaptive | 18/18 | 12 | 4.05 | 1.10 | 0.74 |
| 150 | pre_hybrid | weighted_hint | 18/18 | 12 | 2.87 | 0.00 | 0.73 |
| 150 | pre_hybrid | weighted_nohint | 18/18 | 12 | 2.92 | 0.00 | 0.74 |
| 150 | pre_hybrid | weighted_scipy | 18/18 | 12 | 3.22 | 0.00 | 1.14 |
| 150 | pre_hybrid | weighted_three | 18/18 | 12 | 3.62 | 0.64 | 0.72 |
| 150 | public | equal_nohint | 18/18 | 12 | 3.23 | 0.00 | 1.11 |
| 150 | public | weighted_adaptive | 18/18 | 12 | 6.18 | 2.99 | 1.01 |
| 150 | public | weighted_hint | 18/18 | 12 | 3.33 | 0.00 | 1.16 |
| 150 | public | weighted_nohint | 18/18 | 12 | 3.41 | 0.00 | 1.17 |
| 150 | public | weighted_scipy | 18/18 | 11 | 3.30 | 0.00 | 1.09 |
| 150 | public | weighted_three | 18/18 | 12 | 5.28 | 1.82 | 1.06 |
| 600 | pre_hybrid | equal_nohint | 3/3 | 3 | 433.02 | 0.00 | 425.65 |
| 600 | pre_hybrid | weighted_adaptive | 3/3 | 3 | 437.95 | 33.51 | 397.49 |
| 600 | pre_hybrid | weighted_hint | 3/3 | 3 | 462.09 | 0.00 | 453.96 |
| 600 | pre_hybrid | weighted_nohint | 3/3 | 3 | 439.62 | 0.00 | 432.29 |
| 600 | pre_hybrid | weighted_scipy | 3/3 | 3 | 583.12 | 0.00 | 575.00 |
| 600 | pre_hybrid | weighted_three | 3/3 | 3 | 464.30 | 3.32 | 453.18 |
| 600 | public | equal_nohint | 3/3 | 3 | 461.02 | 0.00 | 452.55 |
| 600 | public | weighted_adaptive | 3/3 | 3 | 482.63 | 57.05 | 417.64 |
| 600 | public | weighted_hint | 3/3 | 3 | 275.46 | 0.00 | 271.71 |
| 600 | public | weighted_nohint | 3/3 | 3 | 465.14 | 0.00 | 457.04 |
| 600 | public | weighted_scipy | 3/3 | 3 | 571.01 | 0.00 | 562.69 |
| 600 | public | weighted_three | 3/3 | 3 | 239.89 | 9.78 | 226.12 |

Medians describe this fixed mixture of fixtures; compare matched individual cases before ranking methods.

# Efficiency ablation

| Case | Start | Variant | Budget | Repeat | C | Lower | Proven | Prep s | LNS s | Global s | Nodes | RSS MiB | Status |
|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| stress_10sp_40chr_12ev_mild | pre_hybrid | equal_nohint | 150 | 0 | 17768 | 4343 | False | 7.09 | 0.00 | 142.65 | 6 | 131.0 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | equal_nohint | 150 | 1 | 17768 | 7943 | False | 3.19 | 0.00 | 146.84 | 6 | 405.4 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | equal_nohint | 600 | 0 | 8377 | 8377 | True | 6.48 | 0.00 | 425.65 | 7 | 423.8 | complete |
| stress_10sp_40chr_12ev_mild | public | equal_nohint | 150 | 0 | 74306 | 4343 | False | 7.08 | 0.00 | 142.70 | 6 | 130.9 | complete |
| stress_10sp_40chr_12ev_mild | public | equal_nohint | 150 | 1 | 74306 | 7925 | False | 3.27 | 0.00 | 146.69 | 6 | 391.5 | complete |
| stress_10sp_40chr_12ev_mild | public | equal_nohint | 600 | 0 | 8377 | 8377 | True | 3.27 | 0.00 | 269.09 | 7 | 417.5 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_adaptive | 150 | 0 | 14092 | 4343 | False | 6.74 | 36.92 | 106.12 | 6 | 151.8 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_adaptive | 150 | 1 | 14092 | 4343 | False | 3.26 | 25.74 | 120.91 | 6 | 151.9 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 6.17 | 33.51 | 397.49 | 7 | 439.5 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_adaptive | 150 | 0 | 17016 | 4343 | False | 7.05 | 36.82 | 105.94 | 6 | 153.4 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_adaptive | 150 | 1 | 16103 | 4343 | False | 3.20 | 35.21 | 111.51 | 6 | 153.7 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 6.55 | 57.56 | 417.64 | 7 | 439.8 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_hint | 150 | 0 | 17768 | 4343 | False | 6.62 | 0.00 | 143.17 | 6 | 131.3 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_hint | 150 | 1 | 17768 | 4343 | False | 3.24 | 0.00 | 146.69 | 6 | 134.2 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_hint | 600 | 0 | 8377 | 8377 | True | 7.18 | 0.00 | 453.96 | 7 | 418.1 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_hint | 150 | 0 | 97002 | 4343 | False | 7.35 | 0.00 | 142.45 | 6 | 133.7 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_hint | 150 | 1 | 74306 | 7854 | False | 3.22 | 0.00 | 146.87 | 6 | 407.6 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_hint | 600 | 0 | 8377 | 8377 | True | 3.28 | 0.00 | 271.71 | 7 | 423.0 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_nohint | 150 | 0 | 17768 | 4343 | False | 7.14 | 0.00 | 142.60 | 6 | 129.6 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_nohint | 150 | 1 | 17768 | 4343 | False | 3.31 | 0.00 | 146.55 | 6 | 134.0 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 6.65 | 0.00 | 454.15 | 7 | 425.1 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_nohint | 150 | 0 | 74306 | 4343 | False | 7.04 | 0.00 | 142.75 | 6 | 133.6 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_nohint | 150 | 1 | 74306 | 7894 | False | 3.21 | 0.00 | 146.86 | 6 | 381.2 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 3.24 | 0.00 | 271.31 | 7 | 421.4 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_scipy | 150 | 0 | 17768 | 4230 | False | 6.80 | 0.00 | 142.91 | 5 | 151.2 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_scipy | 150 | 1 | 17768 | 4343 | False | 3.28 | 0.00 | 146.60 | 6 | 155.9 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 7.06 | 0.00 | 575.00 | 7 | 540.2 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_scipy | 150 | 0 | 97022 | 4230 | False | 7.12 | 0.00 | 142.64 | 5 | 151.2 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_scipy | 150 | 1 | 74306 | 7938 | False | 3.25 | 0.00 | 146.96 | 6 | 378.5 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 7.23 | 0.00 | 562.69 | 7 | 438.8 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_three | 150 | 0 | 17768 | 4343 | False | 7.14 | 3.90 | 138.68 | 6 | 153.1 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_three | 150 | 1 | 17768 | 4343 | False | 3.41 | 1.76 | 144.67 | 6 | 152.8 | complete |
| stress_10sp_40chr_12ev_mild | pre_hybrid | weighted_three | 600 | 0 | 8377 | 8377 | True | 6.77 | 3.32 | 453.18 | 7 | 432.8 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_three | 150 | 0 | 23916 | 4343 | False | 6.78 | 13.12 | 129.92 | 6 | 154.1 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_three | 150 | 1 | 23916 | 7742 | False | 3.22 | 6.55 | 140.15 | 6 | 432.9 | complete |
| stress_10sp_40chr_12ev_mild | public | weighted_three | 600 | 0 | 8377 | 8377 | True | 6.55 | 12.10 | 419.27 | 7 | 440.2 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | equal_nohint | 150 | 0 | 14394 | 4343 | False | 7.13 | 0.00 | 142.38 | 6 | 132.5 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | equal_nohint | 150 | 1 | 14394 | 4343 | False | 3.33 | 0.00 | 146.53 | 6 | 135.4 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | equal_nohint | 600 | 0 | 8377 | 8377 | True | 3.23 | 0.00 | 228.08 | 7 | 427.4 | complete |
| stress_10sp_40chr_12ev_random | public | equal_nohint | 150 | 0 | 74306 | 4343 | False | 7.22 | 0.00 | 142.55 | 6 | 130.2 | complete |
| stress_10sp_40chr_12ev_random | public | equal_nohint | 150 | 1 | 74306 | 4343 | False | 3.29 | 0.00 | 146.59 | 6 | 135.0 | complete |
| stress_10sp_40chr_12ev_random | public | equal_nohint | 600 | 0 | 8377 | 8377 | True | 6.86 | 0.00 | 452.55 | 7 | 417.5 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_adaptive | 150 | 0 | 12346 | 4343 | False | 7.21 | 36.76 | 105.78 | 6 | 152.0 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_adaptive | 150 | 1 | 12346 | 4343 | False | 3.19 | 31.30 | 115.42 | 6 | 153.2 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 3.21 | 24.40 | 211.63 | 7 | 437.5 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_adaptive | 150 | 0 | 13868 | 4343 | False | 7.20 | 36.68 | 105.92 | 6 | 151.7 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_adaptive | 150 | 1 | 13550 | 4343 | False | 3.27 | 37.73 | 108.91 | 6 | 152.3 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 3.24 | 33.02 | 220.90 | 7 | 439.7 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_hint | 150 | 0 | 14394 | 4230 | False | 7.21 | 0.00 | 142.54 | 5 | 131.6 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_hint | 150 | 1 | 14394 | 4343 | False | 3.30 | 0.00 | 146.55 | 6 | 132.5 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_hint | 600 | 0 | 8377 | 8377 | True | 3.21 | 0.00 | 224.19 | 7 | 422.3 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_hint | 150 | 0 | 196253 | 4343 | False | 7.02 | 0.00 | 143.04 | 6 | 166.8 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_hint | 150 | 1 | 277690 | 4343 | False | 3.27 | 0.00 | 146.77 | 6 | 173.7 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_hint | 600 | 0 | 8377 | 8377 | True | 3.30 | 0.00 | 223.53 | 7 | 419.5 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_nohint | 150 | 0 | 14394 | 4230 | False | 6.41 | 0.00 | 143.33 | 5 | 133.1 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_nohint | 150 | 1 | 14394 | 4343 | False | 3.25 | 0.00 | 146.61 | 6 | 135.3 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 3.20 | 0.00 | 227.23 | 7 | 424.0 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_nohint | 150 | 0 | 74306 | 4343 | False | 7.04 | 0.00 | 142.74 | 6 | 133.3 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_nohint | 150 | 1 | 74306 | 4343 | False | 3.25 | 0.00 | 146.65 | 6 | 136.8 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 7.12 | 0.00 | 457.04 | 7 | 413.2 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_scipy | 150 | 0 | 14394 | 4230 | False | 7.40 | 0.00 | 142.28 | 5 | 153.0 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_scipy | 150 | 1 | 14394 | 4343 | False | 3.26 | 0.00 | 146.60 | 6 | 153.9 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 3.21 | 0.00 | 243.59 | 7 | 475.9 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_scipy | 150 | 0 | 180905 | 4230 | False | 6.56 | 0.00 | 143.29 | 5 | 177.1 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_scipy | 150 | 1 | 225514 | 4343 | False | 3.28 | 0.00 | 146.66 | 6 | 176.7 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 3.20 | 0.00 | 245.32 | 7 | 429.8 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_three | 150 | 0 | 14394 | 4343 | False | 7.24 | 3.54 | 139.00 | 6 | 152.0 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_three | 150 | 1 | 14394 | 4343 | False | 3.32 | 1.77 | 144.80 | 6 | 154.0 | complete |
| stress_10sp_40chr_12ev_random | pre_hybrid | weighted_three | 600 | 0 | 8377 | 8377 | True | 3.18 | 1.42 | 211.50 | 7 | 439.3 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_three | 150 | 0 | 19230 | 4343 | False | 6.30 | 16.06 | 127.48 | 6 | 155.6 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_three | 150 | 1 | 19230 | 4343 | False | 3.20 | 10.13 | 136.61 | 6 | 154.3 | complete |
| stress_10sp_40chr_12ev_random | public | weighted_three | 600 | 0 | 8377 | 8377 | True | 3.22 | 8.63 | 222.28 | 7 | 439.8 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | equal_nohint | 150 | 0 | 16605 | 4343 | False | 6.88 | 0.00 | 142.86 | 6 | 130.8 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | equal_nohint | 150 | 1 | 16605 | 4343 | False | 3.34 | 0.00 | 146.52 | 6 | 134.5 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | equal_nohint | 600 | 0 | 8377 | 8377 | True | 7.36 | 0.00 | 494.52 | 7 | 421.9 | complete |
| stress_10sp_40chr_12ev_strong | public | equal_nohint | 150 | 0 | 74306 | 4343 | False | 6.30 | 0.00 | 143.51 | 6 | 131.7 | complete |
| stress_10sp_40chr_12ev_strong | public | equal_nohint | 150 | 1 | 74306 | 4343 | False | 3.24 | 0.00 | 146.66 | 6 | 136.7 | complete |
| stress_10sp_40chr_12ev_strong | public | equal_nohint | 600 | 0 | 8377 | 8377 | True | 7.13 | 0.00 | 490.11 | 7 | 410.3 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_adaptive | 150 | 0 | 11600 | 4343 | False | 7.25 | 36.21 | 106.33 | 6 | 152.6 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_adaptive | 150 | 1 | 11582 | 4343 | False | 3.25 | 33.79 | 112.86 | 6 | 171.3 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 7.12 | 53.16 | 479.17 | 7 | 443.3 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_adaptive | 150 | 0 | 22696 | 4343 | False | 6.95 | 36.77 | 106.10 | 6 | 153.2 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_adaptive | 150 | 1 | 22696 | 4343 | False | 3.34 | 35.08 | 111.48 | 6 | 153.3 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_adaptive | 600 | 0 | 8377 | 8377 | True | 7.19 | 57.05 | 496.53 | 7 | 436.4 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_hint | 150 | 0 | 16605 | 4343 | False | 7.18 | 0.00 | 142.58 | 6 | 131.1 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_hint | 150 | 1 | 16605 | 4343 | False | 3.27 | 0.00 | 146.59 | 6 | 134.2 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_hint | 600 | 0 | 8377 | 8377 | True | 7.15 | 0.00 | 482.79 | 7 | 420.4 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_hint | 150 | 0 | 157471 | 4343 | False | 7.12 | 0.00 | 142.92 | 6 | 173.0 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_hint | 150 | 1 | 238812 | 4343 | False | 3.32 | 0.00 | 146.69 | 6 | 171.5 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_hint | 600 | 0 | 8377 | 8377 | True | 7.16 | 0.00 | 493.99 | 7 | 409.0 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_nohint | 150 | 0 | 16605 | 4343 | False | 7.29 | 0.00 | 142.45 | 6 | 131.3 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_nohint | 150 | 1 | 16605 | 4343 | False | 3.24 | 0.00 | 146.63 | 6 | 134.2 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 6.43 | 0.00 | 432.29 | 7 | 413.2 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_nohint | 150 | 0 | 74326 | 4230 | False | 7.26 | 0.00 | 142.52 | 5 | 131.3 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_nohint | 150 | 1 | 74306 | 4343 | False | 3.27 | 0.00 | 146.60 | 6 | 135.0 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_nohint | 600 | 0 | 8377 | 8377 | True | 7.15 | 0.00 | 489.91 | 7 | 409.4 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_scipy | 150 | 0 | 16605 | 4230 | False | 6.63 | 0.00 | 143.06 | 5 | 151.2 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_scipy | 150 | 1 | 16605 | 4343 | False | 3.27 | 0.00 | 146.59 | 6 | 156.3 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 7.11 | 0.00 | 582.95 | 7 | 481.2 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_scipy | 150 | 0 | 191634 | 4230 | False | 7.39 | 0.00 | 142.49 | 5 | 177.2 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_scipy | 150 | 1 | 194530 | 4343 | False | 3.24 | 0.00 | 146.70 | 6 | 181.7 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_scipy | 600 | 0 | 8377 | 8377 | True | 6.91 | 0.00 | 577.08 | 7 | 433.3 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_three | 150 | 0 | 16605 | 4343 | False | 7.12 | 3.50 | 139.17 | 6 | 152.4 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_three | 150 | 1 | 16605 | 4343 | False | 3.22 | 1.69 | 145.00 | 6 | 152.1 | complete |
| stress_10sp_40chr_12ev_strong | pre_hybrid | weighted_three | 600 | 0 | 8377 | 8377 | True | 7.37 | 3.68 | 488.47 | 7 | 444.2 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_three | 150 | 0 | 25540 | 4343 | False | 7.43 | 18.10 | 124.28 | 6 | 153.0 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_three | 150 | 1 | 25540 | 4343 | False | 3.29 | 11.13 | 135.49 | 6 | 154.0 | complete |
| stress_10sp_40chr_12ev_strong | public | weighted_three | 600 | 0 | 8377 | 8377 | True | 3.37 | 9.78 | 226.12 | 7 | 437.4 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | equal_nohint | 150 | 0 | 341 | 341 | True | 0.21 | 0.00 | 0.20 | 5 | 51.1 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | equal_nohint | 150 | 1 | 341 | 341 | True | 0.19 | 0.00 | 0.21 | 5 | 51.0 | complete |
| stress_6sp_20chr_5ev_mild | public | equal_nohint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.42 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_mild | public | equal_nohint | 150 | 1 | 341 | 341 | True | 0.31 | 0.00 | 0.42 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.21 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.17 | 0.00 | 0.64 | 5 | 50.9 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.34 | 0.00 | 0.43 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.34 | 0.00 | 0.42 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_hint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.21 | 5 | 51.0 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_hint | 150 | 1 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 50.9 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_hint | 150 | 0 | 341 | 341 | True | 0.34 | 0.00 | 0.41 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_hint | 150 | 1 | 341 | 341 | True | 0.34 | 0.00 | 0.42 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 51.5 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.18 | 0.00 | 0.56 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.43 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.46 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_scipy | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.43 | 5 | 70.4 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.20 | 0.00 | 0.46 | 5 | 70.4 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_scipy | 150 | 0 | 341 | 269 | False | 0.35 | 0.00 | 10.55 | 7 | 70.6 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.47 | 8 | 70.6 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_three | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 51.1 | complete |
| stress_6sp_20chr_5ev_mild | pre_hybrid | weighted_three | 150 | 1 | 341 | 341 | True | 0.19 | 0.00 | 0.20 | 5 | 51.5 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_three | 150 | 0 | 341 | 341 | True | 0.35 | 0.00 | 0.44 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_mild | public | weighted_three | 150 | 1 | 341 | 341 | True | 0.34 | 0.00 | 0.43 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | equal_nohint | 150 | 0 | 341 | 341 | True | 0.19 | 0.00 | 1.62 | 5 | 51.5 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | equal_nohint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.55 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_random | public | equal_nohint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.45 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | public | equal_nohint | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 0.61 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.21 | 0.00 | 0.21 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.55 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.34 | 0.00 | 0.45 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 0.61 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_hint | 150 | 0 | 341 | 341 | True | 0.22 | 0.00 | 0.24 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_hint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.56 | 5 | 51.1 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_hint | 150 | 0 | 341 | 341 | True | 0.34 | 0.00 | 0.45 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_hint | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 0.62 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.22 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.56 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.32 | 0.00 | 1.50 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 0.61 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_scipy | 150 | 0 | 341 | 341 | True | 0.21 | 0.00 | 0.43 | 5 | 70.3 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 1.13 | 5 | 70.5 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_scipy | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.45 | 8 | 70.6 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 1.09 | 8 | 70.7 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_three | 150 | 0 | 341 | 341 | True | 0.19 | 0.00 | 0.21 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_random | pre_hybrid | weighted_three | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.56 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_three | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.45 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_random | public | weighted_three | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 0.61 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | equal_nohint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 51.6 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | equal_nohint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.57 | 5 | 51.5 | complete |
| stress_6sp_20chr_5ev_strong | public | equal_nohint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.46 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_strong | public | equal_nohint | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.44 | 8 | 52.1 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.23 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.57 | 5 | 51.3 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_adaptive | 150 | 0 | 341 | 341 | True | 0.34 | 0.00 | 0.44 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_adaptive | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.46 | 8 | 52.4 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_hint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 51.1 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_hint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.56 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_hint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.45 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_hint | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.61 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.21 | 5 | 51.7 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.57 | 5 | 51.4 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_nohint | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.44 | 8 | 52.5 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_nohint | 150 | 1 | 341 | 341 | True | 0.33 | 0.00 | 0.52 | 8 | 52.2 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_scipy | 150 | 0 | 341 | 341 | True | 0.22 | 0.00 | 0.43 | 5 | 70.5 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 1.15 | 5 | 70.3 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_scipy | 150 | 0 | 341 | 341 | True | 0.35 | 0.00 | 1.01 | 8 | 70.7 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_scipy | 150 | 1 | 341 | 341 | True | 0.15 | 0.00 | 1.10 | 8 | 70.7 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_three | 150 | 0 | 341 | 341 | True | 0.20 | 0.00 | 0.20 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_strong | pre_hybrid | weighted_three | 150 | 1 | 341 | 341 | True | 0.09 | 0.00 | 0.57 | 5 | 51.2 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_three | 150 | 0 | 341 | 341 | True | 0.33 | 0.00 | 0.44 | 8 | 52.3 | complete |
| stress_6sp_20chr_5ev_strong | public | weighted_three | 150 | 1 | 341 | 341 | True | 0.34 | 0.00 | 0.46 | 8 | 52.5 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.45 | 0.00 | 0.76 | 8 | 60.9 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | equal_nohint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.38 | 8 | 61.2 | complete |
| stress_8sp_30chr_8ev_mild | public | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.52 | 0.00 | 1.17 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_mild | public | equal_nohint | 150 | 1 | 1258 | 1258 | True | 0.68 | 0.00 | 0.56 | 8 | 61.0 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.51 | 1.08 | 0.72 | 8 | 82.0 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 0.71 | 0.62 | 0.37 | 8 | 81.6 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.55 | 3.14 | 1.11 | 8 | 81.9 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 0.68 | 1.13 | 0.52 | 8 | 82.2 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.51 | 0.00 | 0.73 | 8 | 60.7 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_hint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.37 | 8 | 60.8 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.40 | 0.00 | 1.21 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_hint | 150 | 1 | 1258 | 1258 | True | 0.69 | 0.00 | 0.58 | 8 | 61.4 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.49 | 0.00 | 0.75 | 8 | 61.2 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.38 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.55 | 0.00 | 1.18 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 0.69 | 0.00 | 0.55 | 8 | 61.6 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.40 | 0.00 | 1.62 | 8 | 82.0 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.72 | 8 | 82.2 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.54 | 0.00 | 1.64 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 0.69 | 0.00 | 0.74 | 8 | 82.0 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.49 | 0.66 | 0.72 | 8 | 81.6 | complete |
| stress_8sp_30chr_8ev_mild | pre_hybrid | weighted_three | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.43 | 0.37 | 8 | 81.7 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.50 | 1.81 | 1.14 | 8 | 82.4 | complete |
| stress_8sp_30chr_8ev_mild | public | weighted_three | 150 | 1 | 1258 | 1258 | True | 0.69 | 0.75 | 0.52 | 8 | 82.2 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.60 | 0.00 | 0.75 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | equal_nohint | 150 | 1 | 1258 | 1258 | True | 1.50 | 0.00 | 0.75 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_random | public | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.51 | 0.00 | 1.17 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_random | public | equal_nohint | 150 | 1 | 1258 | 1258 | True | 1.34 | 0.00 | 1.05 | 8 | 61.7 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.48 | 1.17 | 0.77 | 8 | 81.5 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 1.45 | 1.13 | 0.76 | 8 | 81.6 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.51 | 3.01 | 1.09 | 8 | 82.4 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 1.29 | 3.25 | 0.93 | 8 | 81.9 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.41 | 0.00 | 0.73 | 8 | 60.5 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_hint | 150 | 1 | 1258 | 1258 | True | 1.52 | 0.00 | 0.74 | 8 | 60.5 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.55 | 0.00 | 1.19 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_hint | 150 | 1 | 1258 | 1258 | True | 1.51 | 0.00 | 1.14 | 8 | 61.2 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.53 | 0.00 | 0.76 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 1.50 | 0.00 | 0.74 | 8 | 61.3 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.50 | 0.00 | 1.16 | 8 | 61.3 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 1.39 | 0.00 | 1.09 | 8 | 61.6 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.47 | 0.00 | 1.58 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 0.66 | 0.00 | 0.40 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.49 | 0.00 | 1.05 | 8 | 81.9 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 1.49 | 0.00 | 0.89 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.51 | 0.83 | 0.76 | 8 | 81.3 | complete |
| stress_8sp_30chr_8ev_random | pre_hybrid | weighted_three | 150 | 1 | 1258 | 1258 | True | 1.51 | 0.64 | 0.76 | 8 | 81.7 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.45 | 2.03 | 1.06 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_random | public | weighted_three | 150 | 1 | 1258 | 1258 | True | 1.45 | 2.05 | 1.06 | 8 | 82.2 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.52 | 0.00 | 0.74 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | equal_nohint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.37 | 8 | 61.2 | complete |
| stress_8sp_30chr_8ev_strong | public | equal_nohint | 150 | 0 | 1258 | 1258 | True | 1.51 | 0.00 | 1.19 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_strong | public | equal_nohint | 150 | 1 | 1258 | 1258 | True | 0.69 | 0.00 | 0.55 | 8 | 61.7 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.51 | 1.12 | 0.77 | 8 | 81.6 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 1.36 | 0.93 | 0.66 | 8 | 81.9 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_adaptive | 150 | 0 | 1258 | 1258 | True | 1.43 | 2.97 | 1.16 | 8 | 82.4 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_adaptive | 150 | 1 | 1258 | 1258 | True | 0.71 | 1.11 | 0.53 | 8 | 82.3 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.52 | 0.00 | 0.73 | 8 | 60.8 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_hint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.37 | 8 | 60.7 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_hint | 150 | 0 | 1258 | 1258 | True | 1.56 | 0.00 | 1.20 | 8 | 61.6 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_hint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.58 | 8 | 61.1 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.50 | 0.00 | 0.74 | 8 | 61.2 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.38 | 8 | 61.3 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_nohint | 150 | 0 | 1258 | 1258 | True | 1.46 | 0.00 | 1.20 | 8 | 61.5 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_nohint | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.56 | 8 | 61.4 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.46 | 0.00 | 0.93 | 8 | 82.2 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 1.54 | 0.00 | 0.91 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_scipy | 150 | 0 | 1258 | 1258 | True | 1.50 | 0.00 | 0.91 | 8 | 82.1 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_scipy | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.00 | 0.69 | 8 | 82.0 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.46 | 0.64 | 0.73 | 8 | 81.2 | complete |
| stress_8sp_30chr_8ev_strong | pre_hybrid | weighted_three | 150 | 1 | 1258 | 1258 | True | 1.36 | 0.55 | 0.65 | 8 | 81.3 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_three | 150 | 0 | 1258 | 1258 | True | 1.58 | 1.84 | 1.15 | 8 | 82.0 | complete |
| stress_8sp_30chr_8ev_strong | public | weighted_three | 150 | 1 | 1258 | 1258 | True | 0.70 | 0.81 | 0.53 | 8 | 82.5 | complete |

Expected 252 tasks; found 252. Audit case records: 9/9.
Completion is not proof; only matching global bounds count as reported optima.
Equal vs weighted allocation both retain proven-zero components and visit small unresolved models first.
The pre_hybrid start is frozen historical evidence, not the latest known optimum.
Two repeats at 150 s; one per largest case/variant/start at 600 s. Timing differences need repeated confirmation.
Stage timings overlap. Incumbent trajectories are improvements, not global bound trajectories.
AUDIT stress_6sp_20chr_5ev_mild: strict-improvement infeasibility verified
AUDIT stress_6sp_20chr_5ev_strong: strict-improvement infeasibility verified
AUDIT stress_6sp_20chr_5ev_random: strict-improvement infeasibility verified
AUDIT stress_8sp_30chr_8ev_mild: strict-improvement infeasibility verified
AUDIT stress_8sp_30chr_8ev_strong: strict-improvement infeasibility verified
AUDIT stress_8sp_30chr_8ev_random: strict-improvement infeasibility verified
AUDIT stress_10sp_40chr_12ev_mild: re-score passed; proof unresolved within audit budget
AUDIT stress_10sp_40chr_12ev_strong: re-score passed; proof unresolved within audit budget
AUDIT stress_10sp_40chr_12ev_random: re-score passed; proof unresolved within audit budget

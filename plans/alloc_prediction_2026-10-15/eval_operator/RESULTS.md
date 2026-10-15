# Computed results

These are offline predictions, not measured speedups.

## Prediction from 100 draws per chain

Tau ratios use only the following 300 draws and the eight slowest coordinates in those draws. Lag-10 errors use the immediately following 100 draws. No evaluation draw enters the training covariance.

| Dataset | Predictor | Median lag-10 absolute error | Median predicted / future tau, slowest eight |
|---|---|---:|---:|
| PM burn dirs | scalar | 0.146 | 0.139 |
| PM burn dirs | joint_group_0 | 0.083 | 0.672 |
| PM burn dirs | joint_aug_0 | 0.098 | 0.826 |
| PM burn dirs | joint_quad5 | 0.110 | 0.969 |
| pair_new_2 dirs | scalar | 0.156 | 0.148 |
| pair_new_2 dirs | joint_group_0 | 0.087 | 0.376 |
| pair_new_2 dirs | joint_quad5 | 0.089 | 0.413 |
| pair_new_3 dirs | scalar | 0.217 | 0.118 |
| pair_new_3 dirs | joint_group_0 | 0.123 | 0.288 |
| pair_new_3 dirs | joint_quad5 | 0.109 | 0.323 |
| pair_point_2 dirs | scalar | 0.141 | 0.102 |
| pair_point_2 dirs | joint_group_0 | 0.102 | 0.206 |
| pair_point_2 dirs | joint_quad5 | 0.101 | 0.224 |
| pair_point_3 dirs | scalar | 0.144 | 0.234 |
| pair_point_3 dirs | joint_group_0 | 0.114 | 0.474 |
| pair_point_3 dirs | joint_quad5 | 0.105 | 0.517 |
| full_T30 sample | scalar | 0.086 | 0.246 |
| full_T30 sample | joint_group_0 | 0.099 | 0.428 |
| full_T30 sample | joint_aug_0 | 0.107 | 0.459 |
| full_T30 sample | joint_quad5 | 0.109 | 0.581 |
| full_forst sample | scalar | 0.044 | 1.008 |
| full_forst sample | joint_group_0 | 0.044 | 1.089 |
| full_forst sample | joint_aug_0 | 0.044 | 1.108 |
| full_forst sample | joint_quad5 | 0.044 | 1.195 |

## First three PM burn decisions

These use the allocator's mean directions fixed from the checkpoint, with the same basis for every chain.

| Training block | Predictor | Slowest-eight predicted / future-300 tau | Lag-10 absolute error |
|---:|---|---:|---:|
| 1 | scalar | 0.174 | 0.150 |
| 1 | joint_group_0 | 0.402 | 0.086 |
| 1 | joint_aug_0 | 0.491 | 0.089 |
| 1 | joint_quad5 | 0.557 | 0.092 |
| 2 | scalar | 0.097 | 0.178 |
| 2 | joint_group_0 | 0.573 | 0.059 |
| 2 | joint_aug_0 | 0.663 | 0.062 |
| 2 | joint_quad5 | 0.736 | 0.078 |
| 3 | scalar | 0.092 | 0.129 |
| 3 | joint_group_0 | 0.384 | 0.064 |
| 3 | joint_aug_0 | 0.458 | 0.095 |
| 3 | joint_quad5 | 0.576 | 0.117 |

## Decisions from identical saved measurements

Costs are predicted reserved CPU seconds per iteration, averaged over chains. Particle rates are the shared maximum over chains. HMC rates include regime rates.

| State | Variant | Cost | Particle rate | HMC rate | Shear rounds |
|---|---|---:|---:|---:|---:|
| pair_point_2 | baseline | 12.277 | 0.874 | 0.137 | 3.732 |
| pair_point_2 | point | 13.080 | 1.000 | 0.087 | 3.739 |
| pair_point_2 | point_no_credit | 12.631 | 0.941 | 0.061 | 3.744 |
| pair_point_2 | kappa | 12.883 | 1.000 | 0.086 | 3.487 |
| pair_point_2 | kappa_no_credit | 12.464 | 0.941 | 0.060 | 3.491 |
| pair_point_2 | own_subject | 15.176 | 0.897 | 0.838 | 3.744 |
| pair_point_2 | baseline_subject | 12.628 | 0.941 | 0.060 | 3.744 |
| pair_point_2 | conditional | 12.439 | 0.901 | 0.105 | 3.744 |
| pair_point_2 | latest | 12.740 | 0.962 | 0.032 | 3.744 |
| pair_new_2 | baseline | 9.656 | 0.653 | 0.121 | 2.964 |
| pair_new_2 | point | 9.677 | 0.653 | 0.122 | 2.990 |
| pair_new_2 | point_no_credit | 9.995 | 0.653 | 0.128 | 3.439 |
| pair_new_2 | kappa | 9.074 | 0.653 | 0.117 | 2.139 |
| pair_new_2 | kappa_no_credit | 9.085 | 0.653 | 0.117 | 2.155 |
| pair_new_2 | own_subject | 12.933 | 0.650 | 0.881 | 3.439 |
| pair_new_2 | baseline_subject | 9.953 | 0.653 | 0.117 | 3.439 |
| pair_new_2 | conditional | 9.895 | 0.642 | 0.132 | 3.439 |
| pair_new_2 | latest | 10.109 | 0.657 | 0.149 | 3.439 |
| pair_new_3 | baseline | 9.976 | 0.597 | 0.156 | 3.736 |
| pair_new_3 | point | 9.977 | 0.597 | 0.156 | 3.733 |
| pair_new_3 | point_no_credit | 9.975 | 0.597 | 0.156 | 3.731 |
| pair_new_3 | kappa | 9.718 | 0.597 | 0.154 | 3.384 |
| pair_new_3 | kappa_no_credit | 9.718 | 0.597 | 0.154 | 3.384 |
| pair_new_3 | latest | 10.081 | 0.599 | 0.175 | 3.731 |
| pair_point_3 | baseline | 12.950 | 0.843 | 0.158 | 4.298 |
| pair_point_3 | point | 13.275 | 0.907 | 0.085 | 4.258 |
| pair_point_3 | point_no_credit | 13.276 | 0.907 | 0.085 | 4.261 |
| pair_point_3 | kappa | 12.954 | 0.907 | 0.083 | 3.785 |
| pair_point_3 | kappa_no_credit | 12.958 | 0.907 | 0.083 | 3.790 |
| pair_point_3 | own_subject | 16.234 | 0.861 | 0.874 | 4.612 |
| pair_point_3 | baseline_subject | 13.269 | 0.907 | 0.083 | 4.261 |
| pair_point_3 | latest | 13.370 | 0.927 | 0.060 | 4.261 |
| arm_pointsave_1 | baseline | 14.129 | 0.881 | 0.120 | 6.781 |
| arm_pointsave_1 | point | 15.199 | 1.000 | 0.136 | 6.712 |
| arm_pointsave_1 | point_no_credit | 14.462 | 0.921 | 0.111 | 6.732 |
| arm_pointsave_1 | kappa | 13.633 | 1.000 | 0.080 | 4.777 |
| arm_pointsave_1 | kappa_no_credit | 13.055 | 0.921 | 0.053 | 4.777 |
| arm_pointsave_1 | own_subject | 17.438 | 0.911 | 0.908 | 6.737 |
| arm_pointsave_1 | baseline_subject | 14.255 | 0.921 | 0.053 | 6.737 |
| arm_pointsave_1 | conditional | 14.257 | 0.915 | 0.069 | 6.776 |
| arm_pointsave_1 | latest | 14.386 | 0.922 | 0.144 | 6.406 |
| checkpoint | baseline | 10.351 | 0.603 | 0.196 | 5.023 |
| checkpoint | point | 10.652 | 0.611 | 0.203 | 5.366 |
| checkpoint | point_no_credit | 10.374 | 0.611 | 0.203 | 4.907 |
| checkpoint | latest | 10.468 | 0.619 | 0.212 | 4.907 |

# Semi-Blind Channel-and-Signal Estimation for Uplink Massive MIMO With Channel Sparsity

This package contains the official implementation of the **semi-blind channel-and-signal estimation (SCSE) algorithm**  and **simplified SCSE (S-SCSE) algorithm** proposed in the paper: 

> W. Yan and X. Yuan, "Semi-Blind Channel-and-Signal Estimation for Uplink Massive MIMO With Channel Sparsity," in IEEE Access, vol. 7, pp. 95008-95020, 2019, doi: [10.1109/ACCESS.2019.2928092] (https://ieeexplore.ieee.org/document/8766122).

## Introduction

This paper proposed a semi-blind channel-and-signal estimation (SCSE) scheme in which the knowledge of the pilot sequences are integrated into the message passing algorithm for sparse matrix factorization. The SCSE algorithm involves enumeration over all possible user permutations, and so is time-consuming when the number of users is relatively large. To reduce complexity,  the simplified SCSE (S-SCSE) is further proposed to accommodate  systems with a large number of users.

## Code Structure

`main.m`: Set system parameters and output simulation results. 

`System_Model`: Generate system model.

`BiGAMP_EM` : EM algorithm updates parameters.

`BiGAMP_test`: demo BiGAMP algorithm.

`SCSE` and `S-SCSE` :  SCSE algorithm and S-SCSE algorithm respectively. 

`SCSE_result.txt`: Save the simulation results. 

## Citation
```
@article{yan2019semi,
  title={Semi-blind channel-and-signal estimation for uplink massive MIMO with channel sparsity},
  author={Yan, Wenjing and Yuan, Xiaojun},
  journal={IEEE Access},
  volume={7},
  pages={95008--95020},
  year={2019},
  publisher={IEEE}
}
```



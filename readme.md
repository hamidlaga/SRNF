# SRNF
## Numerical Inversion of Square Root Normal Field (SRNF) maps

This code implements the algorithm for the numerical inversion of Square Root Normal Field (SRNF) map, published in           
> [Hamid Laga](https://sites.google.com/view/hamidlaga), Qian Xie, [Ian H Jermyn](https://www.durham.ac.uk/staff/i-h-jermyn/), and [Anuj Srivastava](https://anujsrivastava.com/).
> [**Numerical inversion of SRNF maps for elastic shape analysis of genus-zero surfaces**](https://ieeexplore.ieee.org/abstract/document/7807327).
> IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI) 39(12), pp.2451-2464 (2017).

The program requires Matlab, and has been extensively tested on Windows. It should also run on Mac. 
It uses some codes written in C, which need to be compiled with mex. Pre-compiled files for Windows (32 and 64) and Mac (32 and 64) are included.  

It is standalone and does not require any setup (apart from having Matlab running on Windows). It uses some third-party codes, 
which are included in this distribution with proper license. 

This code was developed for research purposes and is not meant for commercial use. We expect that the code may have some bugs and issues. However, while we will be happy to answer some questions, we may not be able to provide support in trouble-shooting the code. 

## Citation
If you use this code, please make sure you cite: 
```bibtex
@article{laga2017numerical,
  title={Numerical inversion of SRNF maps for elastic shape analysis of genus-zero surfaces},
  author={Laga, Hamid and Xie, Qian and Jermyn, Ian H and Srivastava, Anuj},
  journal={IEEE transactions on pattern analysis and machine intelligence},
  volume={39},
  number={12},
  pages={2451--2464},
  year={2017},
  publisher={IEEE}
}
```
and
```bibtex
@article{jermyn2017elastic, 
title={Elastic shape analysis of three-dimensional objects}, 
author={Jermyn, Ian H and Kurtek, Sebastian and Laga, Hamid and Srivastava, Anuj}, 
journal={Synthesis Lectures on Computer Vision}, 
volume={12}, 
number={1}, 
pages={1--185}, 
year={2017}, 
publisher={Morgan \& Claypool Publishers} 
}
```

## How to use the code:
Scripts are provided to illustrate various usage scenarios. In particular:

test_inverseQ.m
Demonstrates how to invert a SRNF map given an initial surface.


## Licence and Copyright
Please refer to the file LICENSE.

## Contributors
Hamid Laga, Sebastian Kurtek, Qian Xie, Ian H. Jermy, Anuj Srivastava.

## Contact

For any questions, you can contact `hamid.laga@gmail.com`.


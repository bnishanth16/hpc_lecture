# HPSC-Final Report
Name: Nishanth Baskaran

Student ID: 19M15017

Implementation:

1. Navier-Stokes Solver for Cavity flow (C++ version)

  Run the following command in the folder [```cpp_version```](https://github.com/bnishanth16/hpc_lecture/tree/master/final_report/cpp_version)
  ```
  g++ main.cpp && ./a.out && python graph_plot.py
  ```
  ``` cavity_flow.png ``` is saved. 
  
  ![cpp-result](https://github.com/bnishanth16/hpc_lecture/blob/master/final_report/cpp_version/cavity_flow.png)
  
2. Navier-Stokes Solver for Cavity flow (cuda version)

  Run the following command in the folder [```cuda_version```](https://github.com/bnishanth16/hpc_lecture/tree/master/final_report/cuda_version)
  ```
  nvcc main.cu && ./a.out && python graph_plot.py
  ```
  ``` cavity_flow_cuda.png ``` is saved.
  
  ![cuda-result](https://github.com/bnishanth16/hpc_lecture/blob/master/final_report/cuda_version/cavity_flow_cuda.png)

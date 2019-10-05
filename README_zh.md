# CDFCI（坐标下降完全组态相互作用算法）
[![CircleCI](https://circleci.com/gh/quan-tum/CDFCI.svg?style=svg)](https://circleci.com/gh/quan-tum/CDFCI)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/quan-tum/CDFCI)

本软件包是一个简单而高效的CDFCI（Coordinate Descent Full Configuration Interaction, 坐标下降完全组态相互作用算法）算法的实现。本软件包使用C++14标准。

## 算法简介
CDFCI是一个在组态相互作用的框架下的一个计算电子结构基态的高效算法。CDFCI把完全组态相互作用的特征值问题重构成一个等价的无约束非凸优化问题，然后通过自适应的坐标下降方法进行求解。CDFCI也使用了确定的压缩策略。因此，CDFCI可以捕捉到重要的Slater行列式，并且更频繁地更新更重要的坐标。

## 引用
如果您使用了本软件包或者CDFCI算法，请引用以下文章。我们感谢您的支持！

- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J. Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI: 10.1021/acs.jctc.9b00138

在下面的文章里，你可以找到CDFCI的一些变种以及收敛性分析。

- Y. Li, J. Lu and Z. Wang, [Coordinate-wise Descent Methods for Leading Eigenvalue Problem](https://doi.org/10.1137/18M1202505), *SIAM J. Sci. Comput.*, *41(4)*, A2681-A2716 (2019), DOI: 10.1137/18M1202505

## 如何使用
首先克隆本软件仓库。
```
git clone https://github.com/quan-tum/CDFCI.git
cd CDFCI
```
你只需要一个支持C++14的编译器， 比如```g++```或者```icc```就可以编译运行了。我们强烈推荐使用因特尔(Intel&reg;)的编译器```icc```以获得高效率。

### 编译
我们使用```Makefile```来辅助编译。首先，请编辑```Makefile```选择需要的C++编译器（默认调用```g++```)。然后运行```make```命令。
```
make
```
之后会生成两个可执行文件。```cdfci```是单线程版本，```cdfci_omp```是支持OpenMP的多线程版本。

### FCIDUMP
FCIDUMP文件储存了量子系统的轨道之间的积分值。这些积分值通常使用Hartree-Fock计算。它们完全定义了系统的哈密顿量。

很多量子化学的软件包可以生成FCIDUMP文件。比如，[Psi4](http://psicode.org/) 从1.2版本开始提供 [```fcidump``` 函数](http://psicode.org/psi4manual/1.2/api/psi4.driver.fcidump.html)。 关于老版本的 Psi4, [HANDE-QMC](https://github.com/hande-qmc/hande) 提供了一个 [插件](https://github.com/hande-qmc/fcidump)。你可以在 [```test/example/h2o_ccpvdz/README.md```](test/example/h2o_ccpvdz/README.md) 找到我们使用的Psi4输入文件。 [PySCF](https://github.com/pyscf/pyscf) 也提供了 [FCIDUMP dump](https://sunqm.github.io/pyscf/tools.html#module-pyscf.tools.fcidump) 工具。

### 输入文件
我们使用JSON作为输入文件的格式。下面是一个典型的输入文件的例子。
```
{
    "hamiltonian": {
        "fcidump_path": "test/example/h2o_ccpvdz/FCIDUMP",
        "threshold": 1e-13
    },
    "solver":{
        "cdfci": {
            "max_memory": 1,
            "num_iterations": 30000,
            "report_interval": 1000,
            "ref_det_occ": [0, 1, 2, 3, 4, 5, 26, 27, 34, 35],
            "z_threshold": 0,
            "z_threshold_search": false
        }
    }
}
```

- ```hamiltonian``` 需要计算的量子系统的哈密顿量
  - ```fcidump_path```  FCIDUMP文件的路径。必需。
  - ```threshold``` CDFCI会忽略小于FCIDUMP文件中小于此阈值的积分值。默认值为1e-13.
- ```solver```
  - ```cdfci``` CDFCI算法求解器
    - ```max_memory``` 波函数最大可使用的内存（GB）。本软件包会根据此内存值来计算波函数的最大可能大小。注意:此处最大可使用的内存并不一定会被完全使用，因为波函数只支持2的幂次方的大小。建议：使用最大的可用内存以得到最精确的结果。默认值：1。
    - ```num_iterations```  算法迭代步数。CDFCI算法每步更新一个坐标。建议：根据求解系统的大小，从几万到几亿都可能。默认值：30000。
    - ```report_interval``` 程序每隔```report_interval```步输出一次计算的信息。默认值：1000。
    - ```ref_det_occ``` 定义参考Slater行列式占用的轨道。这是CDFCI算法的初值。轨道从0开始计数。默认值：程序从FCIDUMP中自动构造Hartree-Fock态作为参考行列式（初值）。
    - ```z_threshold``` 波函数截断阈值。此阈值越小，那么波函数就会越大，占用内存空间就会越多，同时计算精度就越高。如果值为0，那么程序会收敛到精确值。请参考参考文章中的2.2.2章节。默认值：0。
    - ```z_threshold_search``` 如果设置为true，程序会尝试自动搜索可用的```z_threshold```，以保证波函数可以储存在可用内存中。简而言之，如果波函数过大，程序会把```z_threshold```增大十倍，然后重新运行。警告：此为测试功能，在某些情况下并不保证一定能完美运行。默认值：false。

在```test/example```文件夹下面有更多的例子。

### 如何运行
请使用JSON输入文件的路径作为程序的参数。
```
./cdfci input.json
```
如果输入文件中没有定义参考Slater行列式```ref_det_occ```，程序会通过FCIDUMP文件来尝试确定Hartree-Fock态作为初始值。如果FCIDUMP文件里面有轨道能量（比如Psi4生成的FCIDUMP），那么CDFCI程序会通过填充最低能量的轨道来定义Hartree-Fock态。否则，CDFCI程序会假设轨道已经按照能量从低到高排序。从此初始态开始，CDFCI开始进行迭代直到达到```num_iteration```最大迭代步数为止。

我们建议把```num_iteration```设置尽可能大，当算法收敛到需要的精度的时候终止程序进程。我们建议把```max_memory```设置的尽可能大，以获得尽可能高的精度。CDFCI算法只有一个参数```z_threshold```截断阈值需要调整。我们建议首先设置一个很小的值，比如0。如果波函数过大超过了可用的内存```max_memory```，我们可以提高截断阈值```z_threshold```然后重新运行程序，直到波函数可以被储存。我们也可以设置```z_thershold_search```，让程序自动搜索合适的截断阈值。在绝大多数情况，程序可以自动找到一个合适的阈值。搜索```z_threshold```只会花费很少的时间。

请运行以下命令来运行启用OpenMP的多线程CDFCI程序。
```
export OMP_NUM_THREADS=8  # 设置线程数
./cdfci_omp input.json
```
注意：OpenMP支持仍在测试阶段，程序的效率并没有达到最优。因为支持并行的原因，OpenMP版本的CDFCI的单线程效率低于单线程版本的CDFCI，粗略估计，3线程的OpenMP版本的CDFCI和单线程的CDFCI程序速度接近。

### 输出
CDFCI程序会在标准输出里面输出六列信息。
- ```Iteration``` CDFCI迭代步数。每步更新一个坐标。
- ```Energy``` CDFCI计算的变分能量 (**v<sub>t</sub>**<sup>\*</sup> H **v<sub>t</sub>**) / (**v<sub>t</sub>**<sup>\*</sup> **v<sub>t</sub>**)， 其中t是迭代步数。
- ```dx``` 当前波函数**b**的坐标更新值。
- ```|z|_0``` 波函数**c** = H**b** 中非零元的个数。它对应内存占用量，并且可以被```z_threshold```截断阈值参数来控制。
- ```|H_i|_0``` 波函数**c** = H**b** 中当前更新的坐标个数。如果没有压缩```z_threshold``` = 0, 它就是哈密顿量当前列的非零元的个数。如果```z_threshold``` > 0， 那么它是压缩后保留的非零元的个数。
- ```Time``` CDFCI程序计算时间（秒）。

### 测试
请运行以下命令来编译运行测试集。
```
make test
```
程序会测试8个量子系统。这些例子在```test/example```文件夹里面。测试会用CDFCI来求解这些系统几千步迭代，然后检测求解的能量是否偏离参考值。在[test/README.md](test/README.md)里有更多测试的信息。

## 作者
* [王喆](http://zhewang.pro/)，杜克大学数学系
* [李颖洲](http://yingzhouli.com/)，杜克大学数学系
* [鲁剑锋](https://services.math.duke.edu/~jianfeng/)，杜克大学数学系，物理系，化学系

## 致谢
本项目由美国国家科学基金会DMS-1454939基金和OAC-1450280基金， 以及美国能源部DE-SC0019449基金部分支持。

## 开源协议
本软件包采用BSD 3-Clause开源协议。代码仓库根目录的LICENSE文件为协议文件。
# A Fortran code for PM2.5 2D-Var data assimilation

## 关于
这是一个二维变分同化示例程序，使用FORTRAN语言编写，将PM2.5观测和GEOS-Chem模拟的华北地面PM2.5浓度(ug/m3)进行变分同化。因为只对地面浓度操作，因此是二维变分。

## 文件包括：
1. step01.3dvar.f90 二维变分同化程序。运行数据的输入路径是DATA目录，输出路径是RESULT目录；
2. step02.draw.ncl 绘图程序，读取RESULT的结果绘图并输出到FIGURE目录。此为ncl脚本。
3. readme_cwy.pptx 说明文件
4. run01.sh 一个shell脚本，运行它将编译step01.3dvar.f90并运行程序

## 如何运行：
1. 需要安装 LAPACK软件包，安装方法参见 readme_cwy.pptx。必需
2. 需要安装 NETCDF软件包，安装方法参见 readme_cwy.pptx。若修改了 step01.3dvar.f90程序，不使用 nc 数据格式，可不用此软件包。
3. 运行 ./run01.sh （我的操作环境是Mac终端）

## 说明：
1. 程序没有使用复杂的共轭梯度法求解方程，而是用 lapack软件包的函数求解矩阵。求解的方程参见 readme_cwy.pptx
2. 因计算机内存限制，这个示例程序可能不能处理很大的格点阵。只是一个变分示例程序。
3. FORTRAN程序包含读入观测和模式数据、检查观测点是否在模拟区域内、将模拟值插值到观测点并计算两者差异、构造矩阵方程，求解矩阵。
4. 模式数据是 nc 格式，观测数据是文本格式。可按自己需求修改。修改时需注意FORTRAN程序中将模式线性插值到观测点的部分是否适用于与你的个例。或将你的数据写成我的 nc 格式，注意经纬度顺序，则程序就不用修改。
5. FORTRAN程序中需要注意的参数：
   （1）obserr是观测误差，单位ug/m3；
   （2）moderr是模式误差，单位ug/m3。此值只有在284行为.true.时才起使用，否则使用nc数据中的模式格点标准差（单位也是ug/m3）当背景误差；
   （3）R是correlation length，可随意调节看同化效果；
   （4）inpath1是观测数据路径、离散的站点数据
   （5）inpath2是模拟数据路径、格点化的模式结果，二维背景场
   （6）oupath0是程序输出路径。输出文件obspm25.txt是真正参数同化计算的观测（有些观测落在模拟区域外，不参与同化）；da.nc是格点结果，其中的pm25noda, pm25ysda分别是同化前和同化后的pm2.5浓度(ug/m3)。
   除以上参数外，其它参数无需修改。

## 问题与反馈：
changwy@mail.iap.ac.cn

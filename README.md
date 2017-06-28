pwf - Newton-Raphon 法潮流计算程序（库）
========================================

## 文件清单

./pwf : 潮流计算程序（库）文件夹

./mtx : 矩阵运算文件夹

./minq : 最优潮流计算程序文件夹

./minf : 火电厂最优发电程序（已不用）

./lgpl_2_1.txt : LGPL 许可原文

./License.txt : 本软件的授权许可

./Makefile : 备份管理文件

./README.md : 必读文件（此文件）

## 使用说明

pwf 程序库使用 KLU 进行求解，应首先安装 libklu, libamd, libbtf, libcolamd.

进入 pwf 文件夹，输入 make， 即可编译生成可执行文件 pwf 与程序库 libpwf.a.

进入 minq 文件夹，输入 make，即可生成可执行文件 minq.

## pwf 用法简介

pwf [-agrcnlyj] [-lim] [+tol] [-ooutfile] [infile...]

-a	adjusted load flow   
-g	G-S method   
-r	rectangular N-R method  
-c	write to CDF file  
-n	print the nodal information  
-l	print the line flows  
-y	print the YBus  
-j	print the Jacobian  
-lim	set the iteration limit to lim(INT)  
+tol	set the precision to tol(DOUBLE)  
-ofilename  
	set the output file to filename  
infile...  
	CDF input files  

## pwf 说明

The pwf calculates load flows from infiles, and write the results to outfile.
If no options given, pwf reads from stdin and write to stdout, using the 
default settings: limit 150, tolerance .001, and polar form N-R method.

## 版权声明

pwf, mtx, minq, minf 为自由软件，使用LGPL V2.1授权，详见/License.txt。

本软件中使用的测试数据 ieee14, ieee30, ieee118, ieee300 以及CDF格式说明文件/pwf/res/cdfformat.txt 取自https://www2.ee.washington.edu/research/pstca/，
版权归原作者所有，详见 https://www2.ee.washington.edu/research/pstca/pg_tcaintro.htm

## 作者

李扬
邮箱：liyangpostbox@163.com


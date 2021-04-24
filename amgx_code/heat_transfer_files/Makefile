# Copyright (c) 2011-2016 NVIDIA CORPORATION. All Rights Reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto. Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.

# CUDA Toolkit location
INC_D = -I/usr/local/cuda/include -I/home/divyaprakash/AMGX/base/include
LIB_D = -L/usr/local/cuda/lib64 -L/home/divyaprakash/AMGX/build
LIB_L = -lcusparse -lcublas -lcudart -ldl

# Examples
amgx_capi:
	gcc -O2 -std=c99 heat_transfer.c -c $(INC_D)
	g++ -O2 heat_transfer.o -o heat_transfer $(LIB_D) $(LIB_L) -L. -lamgxsh
# g++ -O2 heat_transfer.o -o heat_transfer $(LIB_D) $(LIB_L) -L../lib -lamgxsh -Wl,-rpath=../lib

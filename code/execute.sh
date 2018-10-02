#!/bin/bash
./jssp_linear abz8 665 1 F > abz8_packing_linear.out
#./jssp_binario abz8 665 1 F > abz8_packing_binario.out
./jssp_linear abz8 665 1 Fe > abz8_machine_linear.out
./jssp_binario abz8 665 1 Fe > abz8_machine_binario.out
./jssp_linear la02 665 1 F > la02_packing_linear.out
./jssp_binario la02 665 1 F > la02_packing_binario.out
./jssp_linear la02 665 1 Fe > la02_machine_linear.out
./jssp_binario la02 665 1 Fe > la02_machine_binario.out
./jssp_linear ta05 1190 1 F > ta05_packing_linear.out
./jssp_binario ta05 1190 1 F > ta05_packing_binario.out
./jssp_linear ta05 1190 1 Fe > ta05_machine_linear.out
./jssp_binario ta05 1190 1 Fe > ta05_machine_binario.out
./jssp_linear ta18 1396 1 F > ta18_packing_linear.out
./jssp_binario ta18 1396 1 F > ta18_packing_binario.out
./jssp_linear ta18 1396 1 Fe > ta18_machine_linear.out
./jssp_binario ta18 1396 1 Fe > ta18_machine_binario.out
./jssp_linear dmu01 2563 1 F > dmu01_packing_linear.out
./jssp_binario dmu01 2563 1 F > dmu01_packing_binario.out
./jssp_linear dmu01 2563 1 Fe > dmu01_machine_linear.out
./jssp_binario dmu01 2563 1 Fe > dmu01_machine_binario.out


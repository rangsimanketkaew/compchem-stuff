# CP2K on Piz Daint

Tested version: 9.1

## Install toolchains
```sh
cd cp2k-9.1/tools/toolchain/
./install_cp2k_toolchain.sh -j 16 --with-gcc=system --with-plumed
```

## Compile CP2K

```sh
cd cp2k-9.1
source arch/CRAY-XC40-gnu.psmp 10.3.0
make -j 16 ARCH=CRAY-XC40-gnu VERSION=psmp
```

Executables will be created in `exe/` folder

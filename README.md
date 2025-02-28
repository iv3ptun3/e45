k18geant4
=========

K1.8 geant4 simulation tool.



## Platform

This tool is developed on the platform of local in MacOS(arm64)
- clang++
- ROOT 6.32.08
- Geant4 11.2.1



## How to install

Environment variables should be set.

```shell
export G4WORKDIR=$HOME/work
export PATH=$G4WORKDIR/bin/Darwin-clang:$PATH
. /sw/packages/root/6.32.08/bin/thisroot.sh
. /sw/packages/geant4/11.2.1/bin/geant4.sh
. /sw/packages/geant4/11.2.1/share/Geant4-10.4.3/geant4make/geant4make.sh
```

then

```shell
git clone https://github.com/iv3ptun3/e42.git
cd e42
git checkout main
make
```

for update

```shell
git pull origin main
```

## How to use

Arguments of ConfFile and OutputName are necessary.
G4Macro is an optional argument.

```shell
hyptpc1 [ConfFile] [OutputName] (G4Macro)
hyptpc1 param/conf/default.conf test.root
hyptpc1 param/conf/default.conf test.root g4macro/run.mac
```



## Parameters

Some parameter files that are out of the git control should be linked.
Or scp from kek.cc

```shell
ln -s /group/had/sks/E42/software/param/BEAM/* param/BEAM/
ln -s /group/had/sks/E42/software/param/JAM/* param/JAM/
ln -s /group/had/sks/E42/software/fieldmap .
```


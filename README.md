k18geant4
=========

K1.8 geant4 simulation tool.



## Platform

This tool is developed on the platform of KEKCC, Scientific Linux 6.10.
- g++ (GCC) 4.8.5
- ROOT 6.16.00
- Geant4 10.5.0-mt



## How to install

Environment variables should be set.

```shell
module load gcc/485
module load git/2171
module load python/2.7
module load python/3.5
export G4WORKDIR=$HOME/work/geant4
export PATH=$G4WORKDIR/bin/Linux-g++:$PATH
. /sw/packages/root/6.16.00/bin/thisroot.sh
. /sw/packages/geant4/10.5.0-mt/bin/geant4.sh
. /sw/packages/geant4/10.5.0-mt/share/Geant4-10.5.0/geant4make/geant4make.sh
```

then

```shell
git clone ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18geant4.git
cd k18geant4
git checkout e42
make
```



## How to use

Arguments of ConfFile and OutputName are necessary.
G4Macro is an optional argument.

```shell
hyptpc1 [ConfFile] [OutputName] (G4Macro)
hyptpc1 param/conf/default.conf hoge.root
hyptpc1 param/conf/default.conf hoge.root g4macro/run.mac
```



## Parameters

Some parameter files that are out of the git control should be linked.

```shell
ln -s /group/had/sks/E42/software/param/BEAM/* param/BEAM/
ln -s /group/had/sks/E42/software/param/JAM/* param/JAM/
ln -s /group/had/sks/E42/software/fieldmap .
```


# git clone git@github.com:ablab/algorithmic-biology.git
rm -r build
rm spades_1.0.0*
cp -r assembler build

cd build
cp makefileDebian makefile
cat makefile
debuild -us -uc

cd ..
scp spades_1.0.0* 192.168.222.223:


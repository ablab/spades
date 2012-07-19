mkdir -p bin
cd bin
rm -rf src
rm -rf ext
cp -r ../src .
cp -r ../ext .
cd src
cmake -G "Unix Makefiles" . .
make hammer
cp hammer/hammer ..
make spades
cp debruijn/spades ..
cd ..
rm -r ext
rm -r src


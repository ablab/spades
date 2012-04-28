rm -rf release$(cat assembler/VERSION)
rm -rf build

mkdir -p release$(cat assembler/VERSION)/bayeshammer/

cp -r assembler build
cd build
src/tools/spades_pipeline/prebuild_spades.py 3 199

cd ..

cp -r build/prebuild_spades/bin release$(cat assembler/VERSION)

./make-spades-debian.sh


cp spades_* release$(cat assembler/VERSION)
cp spades.spec release$(cat assembler/VERSION)


scp -r release$(cat assembler/VERSION) builder@192.168.222.254:
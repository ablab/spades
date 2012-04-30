set -e -x
rm -r rpmbuild
rpmdev-setuptree
scp builder@192.168.222.254:release$1/spades.spec .
cp spades.spec rpmbuild/SPECS
cd rpmbuild/SOURCES   
scp builder@192.168.222.254:release$1/spades_$1.tar.gz .
tar -xvf spades_$1.tar.gz
cd spades-$1
touch configure
chmod a+x configure 
cd ..
rm spades_*.tar.gz
tar -czvf spades_$1.tar.gz spades-$1
cd ../SPECS
rpmbuild --define 'version '$1 -ba spades.spec
scp ../RPMS/x86_64/spades-$1-1.fc16.x86_64.rpm builder@192.168.222.254:rpm
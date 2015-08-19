echo $#
echo asd$1asd
timestamp="`date +%Y-%m-%d_%H-%M-%S`"
if [ "$#" > 0 ]
then
    timestamp=$timestamp"_"$1
fi
echo $timestamp 

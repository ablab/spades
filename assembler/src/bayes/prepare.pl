
my $pref = shift;

print $pref . "_scaffolds.fa\n";

system "mkdir $pref";
system "grep -v '^[[:space:]]*>' " . $pref . "_scaffolds.fa > $pref/nocomm.fa";
system "dd if=$pref/nocomm.fa of=$pref/upper.fa conv=ucase";
system "sed 's/N//g' $pref/upper.fa > $pref/noN.fa";
system "cat name.txt $pref/noN.fa > $pref/$pref.fa";
system "rm -rf $pref/noN.fa $pref/nocomm.fa $pref/upper.fa";
system "/home/snikolenko/algorithmic-biology/assembler/src/libs/bowtie-0.12.7/bowtie-build $pref/$pref.fa $pref/$pref";

die "OK.";

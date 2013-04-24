#all 18 species is 0
cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0; foreach $f(1..$#fi){   if ($fi[$f]>0){$tot++;}} if ($tot == 0){print "yes\n";}  ' | wc -l

#primate specific (in all the 4 primate species but not elsewhere)
cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0; $spy_human = 0; $spy_chimp = 0; $spy_macaque = 0; $spy_marmoset = 0; $spy_orangutan = 0; foreach $f(1..$#fi){if ($f==8){if($fi[$f]>0){$spy_human++;}next;}   if ($f==2){if($fi[$f]>0){$spy_chimp++;}next;} if ($f==9){if($fi[$f]>0){$spy_macaque++;}next;} if ($f==10){if($fi[$f]>0){$spy_marmoset++;}next;} if ($f==13){if($fi[$f]>0){$spy_orangutan++;}next;}   if ($fi[$f]>0){$tot++;}} if (($tot == 0) and ($spy_human == 1) and ($spy_chimp == 1) and ($spy_macaque == 1) and ($spy_marmoset == 1) and ($spy_orangutan == 1)){print "$_\n";} '

#primate specific "at least 1" (in at least one of the other primate species but not elsewhere)
cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0; $spy_human = 0; $spy_chimp = 0; $spy_macaque = 0; $spy_marmoset = 0; $spy_orangutan = 0; foreach $f(1..$#fi){if ($f==8){if($fi[$f]>0){$spy_human++;}next;}   if ($f==2){if($fi[$f]>0){$spy_chimp++;}next;} if ($f==9){if($fi[$f]>0){$spy_macaque++;}next;} if ($f==10){if($fi[$f]>0){$spy_marmoset++;}next;} if ($f==13){if($fi[$f]>0){$spy_orangutan++;}next;}   if ($fi[$f]>0){$tot++;}} if (($tot == 0) and ($spy_human == 1) and (($spy_chimp == 1) or ($spy_macaque == 1) or ($spy_marmoset == 1) or ($spy_orangutan == 1))){print "$_\n";} ' | cut -f 1 -d "," > primateSpecific_atLeast1.id

#lncRNA not found in human but found somewhere else (many lower case)
cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0;  $spy_human = 0;  foreach $f(1..$#fi){if ($f==8){if($fi[$f]>0){$spy_human++;}next;} if ($fi[$f]>0){$tot++;}} if (($tot > 0) and ($spy_human == 0) ){print "$_\n";}  '

#lncRNA found in human and in at least one other species
cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0;  $spy_human = 0;  foreach $f(1..$#fi){if ($f==8){if($fi[$f]>0){$spy_human++;}next;} if ($fi[$f]>0){$tot++;}} if (($tot > 0) and ($spy_human == 1) ){print "$_\n";}  '  |wc -l

#lncRNA found in human but not somewhere else
 cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0;  $spy_human = 0;  foreach $f(1..$#fi){if ($f==8){if($fi[$f]>0){$spy_human++;}next;} if ($fi[$f]>0){$tot++;}} if (($tot == 0) and ($spy_human == 1) ){print "$_\n";}  '  |wc -l

#To have the number of lncRNA conserved in a certain number of species you can do this loop
for X in {0..18}; do cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0; foreach $f(1..$#fi){   if ($fi[$f]>0){$tot++;}} if ($tot == '$X'){print "yes\n";}  ' | wc -l ; done
#if you wanna make it cumulative you just have to change == with >= like this:
 for X in {0..18}; do cat heatMap.csv | perl -ne 'chomp $_; if($_=~/chimp/){next;} @fi=split(/,/,$_); $tot=0; foreach $f(1..$#fi){   if ($fi[$f]>0){$tot++;}} if ($tot >= '$X'){print "yes\n";}  ' | wc -l ; done

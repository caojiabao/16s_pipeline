#taxa_census_new.pl

use strict;
use warnings;
use strict;
use Cwd;

my %blacklist=
(
"Clostridiales"=>"unclassified_Clostridiales", 
"Enterobacteriaceae"=>"unclassified_Enterobacteriaceae",
"Erysipelotrichaceae"=>"unclassified_Erysipelotrichaceae",
"Erysipelotrichaceae_incertae_sedis"=>"unclassified_Erysipelotrichaceae", 
"Gp1"=>"unclassified_Gp1",
"Lachnospiraceae"=>"unclassified_Lachnospiraceae", 
"Lachnospiracea_incertae_sedis"=>"unclassified_Lachnospiraceae",
"Peptostreptococcaceae"=>"unclassified_Peptostreptococcaceae", 
"Rhizomicrobium"=>"unclassified_Rhizomicrobium",
"Ruminococcaceae"=>"unclassified_Ruminococcaceae",
"TM7"=>"unclassified_TM7",
"TM7_genera_incertae_sedis"=>"unclassified_TM7",
"OD1"=>"unclassified_OD1", 
"OD1_genera_incertae_sedis"=>"unclassified_OD1"
);


    
    
    my %groups; 
    my @tax; 
    
    my $s=$ARGV[0]; 
    {
    chomp $s; 
    
    {
    open (H, $s);
    my @a=<H>;
    close H; 
    foreach my $line (@a)
    {
    chomp $line;
    push @tax, $line; 

    
    my @ele=split /\t/, $line; 
    my @s=split /\_SEQ/, $ele[0]; 
    if($groups{$s[0]})
    {push @{$groups{$s[0]}},$ele[0]; }
    else{my @temp;
    push @temp, $ele[0];
    $groups{$s[0]}=\@temp; }
    }    
    }
    }
    
    
    
    
    my %first;
    my %second;
    my %third;
    my %fourth;
    my %fifth;
    my %sixth;
    my %seventh;
    my %eighth;
    my %first_otu;
    my %second_otu;
    my %third_otu;
    my %fourth_otu;
    my %fifth_otu;
    my %sixth_otu;
    my %seventh_otu;
    my %eighth_otu;
    #go on... 
    
    
    my @id= sort keys %groups; 
    
    
    
    #for (my $i=0; $tax[$i]; $i+=1)
    while (my $taxi=shift @tax)
    {
    chomp $taxi;
    my @ele=split /\t/, $taxi; 
    $taxi=~ s/\"//gi; 
 $taxi=~ s/ /_/gi;
    $taxi=~ s/\t\t/\t/gi; 
    $taxi=~ s/\///g;
    $taxi=~ /\t(\w+)\tdomain/;
    my $domain=$1; 
    $taxi=~ /\t(\w+)\tphylum/;
    my $phylum=$1;
    $taxi=~ /\t(\w+)\tclass/;
    my $class=$1;
    $taxi=~ /\t(\w+)\torder/;
    my $order=$1;
    #$taxi=~ /\t(\w+)\tfamily/;
$taxi=~ /order\t\d.\d+\t(.+)\tfamily/;    
my $family=$1;
    $taxi=~ /family\t\d.\d+\t(.+)\tgenus/;
    my $genus=$1;
    $taxi=~ /genus\t(\d.\d+)/;
    my $confi=$1; 
    #print $domain, "\t", $phylum, "\t", $class, "\t", $order, "\t", $family, "\t", $genus, "\t", $confi, "\n"; 
    
    if ($domain)
    {
    $first{$ele[0]}=$domain;
    $first_otu{$domain}=0;}
    if ($phylum)
    {
    $second{$ele[0]}=$phylum;
    $second_otu{$phylum}=0;}
    if ($class)
    {
    $third{$ele[0]}=$class;
    $third_otu{$class}=0;}
    if ($order)
    {
    $fourth{$ele[0]}=$order;
    $fourth_otu{$order}=0;}
    if ($family)
    {
    $fifth{$ele[0]}=$family;
    $fifth_otu{$family}=0;}
    
    if ($genus)
    {
    if ($confi>0.8)
    {
    if ($blacklist{$genus})
    {
     
     $sixth{$ele[0]}=$blacklist{$genus};
    $sixth_otu{$blacklist{$genus}}=0;
    }
    else { $sixth{$ele[0]}=$genus;
    $sixth_otu{$genus}=0;}
    }
    else
    {
    my @temp=split /\_/, $family; 
    my $un="unclassified_".$temp[0]; 
    
    $sixth{$ele[0]}=$un;
    $sixth_otu{$un}=0;
    }
    }
    }
    
    
    #die; 
    #ok, create matrix now.... 
    #here comes the boring part......
    #first level
    #----------
    if (%first_otu)
    {
    my @matrix_first1; 
    $matrix_first1[0]="\t";
    my @matrix_first2;
    $matrix_first2[0]="\t";
    my $m=1;
    my $l=1;
    foreach my $otu (sort keys %first_otu)
    {
    $matrix_first1[$m]=$otu."\t";
    $m+=1;
    $matrix_first2[0].=$otu."\t"; 
    }
    
    foreach my $group (@id)
    {
    $matrix_first1[0].=$group."\t";
    $matrix_first2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %first_otu)
    {
    $otu_count{$otu}=0;
    }
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $first{$sequence}, "\n";
    $otu_count{$first{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_first1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_first2[$l].=$otu_count{$otu2}."\t";
    }
    
    $l+=1;
    }
    foreach my $row (@matrix_first1)
    {$row.="\n"};
    foreach my $row2 (@matrix_first2)
    {$row2.="\n"};
    open (FILEHANDLE3, ">kingdom.otu_matrix_column");
    print FILEHANDLE3 @matrix_first1; 
    close FILEHANDLE3;
    
    open (FILEHANDLE4, ">kingdom.otu_matrix_row");
    print FILEHANDLE4 @matrix_first2; 
    close FILEHANDLE4;
    }
    #-----------------------
    
    if (%second_otu)
    {
    my @matrix_second1; 
    $matrix_second1[0]="\t";
    my @matrix_second2;
    $matrix_second2[0]="\t";
    my $m=1;
    foreach my $otu (sort keys %second_otu)
    {
    $matrix_second1[$m]=$otu."\t";
    $m+=1;
    $matrix_second2[0].=$otu."\t"; 
    }
    my $l=1;
    foreach my $group (@id)
    {
    $matrix_second1[0].=$group."\t";
    $matrix_second2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %second_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $second{$sequence}, "\n";
    $otu_count{$second{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_second1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_second2[$l].=$otu_count{$otu2}."\t";
    }
    $l+=1;
    }
    foreach my $row (@matrix_second1)
    {$row.="\n"};
    foreach my $row2 (@matrix_second2)
    {$row2.="\n"};
    
    open (FILEHANDLE5, ">phyla.otu_matrix_column");
    print FILEHANDLE5 @matrix_second1; 
    close FILEHANDLE5;
    
    open (FILEHANDLE6, ">phyla.otu_matrix_row");
    print FILEHANDLE6 @matrix_second2; 
    close FILEHANDLE6;
    }
    #---------------
    
    if (%third_otu)
    {
    my @matrix_third1; 
    $matrix_third1[0]="\t";
    my @matrix_third2;
    $matrix_third2[0]="\t";
    my $m=1;
    foreach my $otu (sort keys %third_otu)
    {
    $matrix_third1[$m]=$otu."\t";
    $m+=1;
    $matrix_third2[0].=$otu."\t"; 
    }
    my $l=1;
    foreach my $group (@id)
    {
    $matrix_third1[0].=$group."\t";
    $matrix_third2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %third_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $third{$sequence}, "\n";
    $otu_count{$third{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_third1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_third2[$l].=$otu_count{$otu2}."\t";
    }
    $l+=1;
    }
    foreach my $row (@matrix_third1)
    {$row.="\n"};
    foreach my $row2 (@matrix_third2)
    {$row2.="\n"};
    
    open (FILEHANDLE7, ">class.otu_matrix_column");
    print FILEHANDLE7 @matrix_third1; 
    close FILEHANDLE7;
    
    open (FILEHANDLE8, ">class.otu_matrix_row");
    print FILEHANDLE8 @matrix_third2; 
    close FILEHANDLE8;
    }
    #--------
    
    if (%fourth_otu)
    {
    my @matrix_fourth1; 
    $matrix_fourth1[0]="\t";
    my @matrix_fourth2;
    $matrix_fourth2[0]="\t";
    my $m=1;
    foreach my $otu (sort keys %fourth_otu)
    {
    $matrix_fourth1[$m]=$otu."\t";
    $m+=1;
    $matrix_fourth2[0].=$otu."\t"; 
    }
    my $l=1;
    foreach my $group (@id)
    {
    $matrix_fourth1[0].=$group."\t";
    $matrix_fourth2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %fourth_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $fourth{$sequence}, "\n";
    $otu_count{$fourth{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_fourth1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_fourth2[$l].=$otu_count{$otu2}."\t";
    }
    $l+=1;
    }
    foreach my $row (@matrix_fourth1)
    {$row.="\n"};
    foreach my $row2 (@matrix_fourth2)
    {$row2.="\n"};
    
    open (FILEHANDLE9, ">order.otu_matrix_column");
    print FILEHANDLE9 @matrix_fourth1; 
    close FILEHANDLE9;
    
    open (FILEHANDLE10, ">order.otu_matrix_row");
    print FILEHANDLE10 @matrix_fourth2; 
    close FILEHANDLE10;
    }
    
    #----------------
    if (%fifth_otu)
    {
    my @matrix_fifth1; 
    $matrix_fifth1[0]="\t";
    my @matrix_fifth2;
    $matrix_fifth2[0]="\t";
    my $m=1;
    foreach my $otu (sort keys %fifth_otu)
    {
    $matrix_fifth1[$m]=$otu."\t";
    $m+=1;
    $matrix_fifth2[0].=$otu."\t"; 
    }
    my $l=1;
    foreach my $group (@id)
    {
    $matrix_fifth1[0].=$group."\t";
    $matrix_fifth2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %fifth_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $fifth{$sequence}, "\n";
    $otu_count{$fifth{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_fifth1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_fifth2[$l].=$otu_count{$otu2}."\t";
    }
    $l+=1;
    }
    foreach my $row (@matrix_fifth1)
    {$row.="\n"};
    foreach my $row2 (@matrix_fifth2)
    {$row2.="\n"};
    open (FILEHANDLE11, ">family.otu_matrix_column");
    print FILEHANDLE11 @matrix_fifth1; 
    close FILEHANDLE11;
    open (FILEHANDLE12, ">family.otu_matrix_row");
    print FILEHANDLE12 @matrix_fifth2; 
    close FILEHANDLE12;
    }
    #---------------------
    if (%sixth_otu)
    {
    my @matrix_sixth1; 
    $matrix_sixth1[0]="\t";
    my @matrix_sixth2;
    $matrix_sixth2[0]="\t";
    my $m=1;
    foreach my $otu (sort keys %sixth_otu)
    {
    $matrix_sixth1[$m]=$otu."\t";
    $m+=1;
    $matrix_sixth2[0].=$otu."\t"; 
    }
    my $l=1;
    foreach my $group (@id)
    {
    $matrix_sixth1[0].=$group."\t";
    $matrix_sixth2[$l]=$group."\t";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %sixth_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $sixth{$sequence}, "\n";
    $otu_count{$sixth{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_sixth1[$i].=$otu_count{$otu2}."\t";
    $i+=1;
    $matrix_sixth2[$l].=$otu_count{$otu2}."\t";
    }
    $l+=1;
    }
    foreach my $row (@matrix_sixth1)
    {$row.="\n"};
    foreach my $row2 (@matrix_sixth2)
    {$row2.="\n"};
    
    
    open (FILEHANDLE13, ">genera.otu_matrix_column");
    print FILEHANDLE13 @matrix_sixth1; 
    close FILEHANDLE13;
    
    open (FILEHANDLE14, ">genera.otu_matrix_row");
    print FILEHANDLE14 @matrix_sixth2; 
    close FILEHANDLE14;
    }
    #--------------
    if (%seventh_otu)
    {
    my @matrix_seventh1; 
    $matrix_seventh1[0]="  ";
    my @matrix_seventh2;
    $matrix_seventh2[0]="  ";
    my $m=1;
    foreach my $otu (sort keys %seventh_otu)
    {
    $matrix_seventh1[$m]=$otu." ";
    $m+=1;
    $matrix_seventh2[0].=$otu." "; 
    }
    my $l=1;
    foreach my $group (keys %groups)
    {
    $matrix_seventh1[0].=$group." ";
    $matrix_seventh2[$l]=$group." ";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %seventh_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $seventh{$sequence}, "\n";
    $otu_count{$seventh{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_seventh1[$i].=$otu_count{$otu2}." ";
    $i+=1;
    $matrix_seventh2[$l].=$otu_count{$otu2}." ";
    }
    $l+=1;
    }
    foreach my $row (@matrix_seventh1)
    {$row.="\n"};
    foreach my $row2 (@matrix_seventh2)
    {$row2.="\n"};
    
    open (FILEHANDLE15, ">seventh.otu_matrix_column");
    print FILEHANDLE15 @matrix_seventh1; 
    close FILEHANDLE15;
    
    open (FILEHANDLE16, ">seventh.otu_matrix_row");
    print FILEHANDLE16 @matrix_seventh2; 
    close FILEHANDLE16;
    }
    #-------------
    if (%eighth_otu)
    {
    my @matrix_eighth1; 
    $matrix_eighth1[0]="  ";
    my @matrix_eighth2;
    $matrix_eighth2[0]="  ";
    my $m=1;
    foreach my $otu (sort keys %eighth_otu)
    {
    $matrix_eighth1[$m]=$otu." ";
    $m+=1;
    $matrix_eighth2[0].=$otu." "; 
    }
    my $l=1;
    foreach my $group (keys %groups)
    {
    $matrix_eighth1[0].=$group." ";
    $matrix_eighth2[$l]=$group." ";
    print $group."\n"; 
    my %otu_count;
    foreach my $otu (keys %eighth_otu)
    {
    $otu_count{$otu}=0;}
    foreach my $sequence (all_seq ($group))
    {
    #print $sequence, "\n";
    ##print "the otu is ", $eighth{$sequence}, "\n";
    $otu_count{$eighth{$sequence}}+=1;
    }
    
    my $i=1;
    foreach my $otu2 (sort keys %otu_count)
    {
    $matrix_eighth1[$i].=$otu_count{$otu2}." ";
    $i+=1;
    $matrix_eighth2[$l].=$otu_count{$otu2}." ";
    }
    $l+=1;
    }
    foreach my $row (@matrix_eighth1)
    {$row.="\n"};
    foreach my $row2 (@matrix_eighth2)
    {$row2.="\n"};
    
    open (FILEHANDLE17, ">eighth.otu_matrix_column");
    print FILEHANDLE17 @matrix_eighth1; 
    close FILEHANDLE17;
    
    open (FILEHANDLE18, ">eighth.otu_matrix_row");
    print FILEHANDLE18 @matrix_eighth2; 
    close FILEHANDLE18;
    }
    
    
    
    #-----------------
    
    
    sub all_seq{
    
    my ($groupname)=@_;
    
    #chomp $group_name;
    #print $group_name,"\n";
    
    my @temp;
    foreach my $t (@{$groups{$groupname}})
    {
        if($first{$t}){push @temp, $t; }
        }
    
        return @temp;
        }
    
    
  
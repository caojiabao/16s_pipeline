
use strict;
use List::Util qw(shuffle);

  my $file = $ARGV[0];

            open(LIST,$file);

            my %hash;
            my $element;
            while (my $line = <LIST>){
                  if ($line =~ />.*/) {

                              $element=$line;
                              $hash{$element}="";
                                    }
                                          else {
                              $hash{$element}.=$line;
                                                      }
                                                      }
                                                      close(LIST);

                                                      my $LINES=$ARGV[1];
                                                      my @ele=split /\./, $file;
                                                      my @out=shuffle (values %hash);
                                                      my $num=@out; 
                                                      for (my $i=1; $i<=$LINES && $i<$num; $i+=1) {
                                                      my $id=$ele[0]."_SEQ".$i;
                                                                                                            print ">".$id."\n".$out[$i];}
                                                       
#   my $index=$#peptides;
                                                        #      push (@out, $peptides[int(rand($index))]);
                                                         #         splice(@peptides, $index, 1);
                                                         #  }

                                                                  #open (OUT, '>out.fasta');
                                                         

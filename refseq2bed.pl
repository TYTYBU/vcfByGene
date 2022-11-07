#!/usr/bin/perl

$infile = $ARGV[0]; # refseq all fields txt file
$infile2 = $ARGV[1]; # pvcf_blocks txt file
$out_path = $ARGV[2]; # output BED file path

if ($out_path eq ""){
  $out_path = ".";
}

open IN, "<", $infile or die;
open IN2, "<", $infile2 or die;

while (my $line = <IN2>){
  chomp $line;
  my ($ind, $chr, $blk, $start_pos, $end_pos) = split /\t/, $line;
  if ($chr == 23){
    $chr = "cX";
  } elsif ($chr == 24){
    $chr = "cY";
  } else {
    $chr = "c" . $chr;
  }
  $blk = "b" . $blk;

  $hash{$chr}{$blk}{start} = $start_pos;
  $hash{$chr}{$blk}{end} = $end_pos;
}

<IN>;
while (my $line = <IN>){
  if($line =~ /[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)/){
    my $gene_id = $1;
    my $chr = $2;
    my $txStart = $3;
    my $txEnd = $4;
    my $exonStart_str = $5;
    my $exonEnd_str = $6;
    my $gene_symbol = $7;

    my $chr2 = substr $chr, 3;
    $chr2 = "c" . $chr2;

    foreach my $blk (keys %{$hash{$chr2}}){
      if ((($hash{$chr2}{$blk}{start} < $txEnd) and ($hash{$chr2}{$blk}{end} > $txEnd)) or (($hash{$chr2}{$blk}{start} < $txStart) and ($hash{$chr2}{$blk}{end} > $txStart))){
        $exonStart_str = substr $exonStart_str, 0, -1;
        @exonStart = split /,/, $exonStart_str;
        $exonEnd_str = substr $exonEnd_str, 0, -1;
        @exonEnd = split /,/, $exonEnd_str;

        open OUT, ">", $out_path . "/" . $chr2 . "_" . $blk . "_" . $gene_symbol . "_" . $gene_id . ".bed";
        foreach my $i (0 .. (scalar(@exonStart)-1)){
          # 5nt flaking regions to include spicing sites for each exon
          print OUT $chr . "\t" . $exonStart[$i]-5 . "\t" . $exonEnd[$i]+5 . "\t" . $gene_symbol . "\n";
        }
        close OUT;

        last;
      }
    }
  }
}

#!/usr/bin/perl

$infile = $ARGV[0]; # mane all fields txt file
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
  if($line =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)/){
    my $chr = $1;
    my $gene_start = $2;
    my $gene_end = $3;
    my $block_size_str = $4;
    my $block_start_str = $5;
    my $gene_symbol = $6;

    my $chr2 = substr $chr, 3;
    $chr2 = "c" . $chr2;

    foreach my $blk (keys %{$hash{$chr2}}){
      if ((($hash{$chr2}{$blk}{start} < $gene_end) and ($hash{$chr2}{$blk}{end} > $gene_end)) or (($hash{$chr2}{$blk}{start} < $gene_start) and ($hash{$chr2}{$blk}{end} > $gene_start))){
        $block_size_str = substr $block_size_str, 0, -1;
        @block_size = split /,/, $block_size_str;
        $block_start_str = substr $block_start_str, 0, -1;
        @block_start = split /,/, $block_start_str;

        open OUT, ">", $out_path . "/" . $chr2 . "_" . $blk . "_" . $gene_symbol . ".bed";
        foreach my $i (0 .. (scalar(@block_size)-1)){
          print OUT $chr . "\t" . ($gene_start+$block_start[$i]) . "\t" . ($gene_start+$block_start[$i]+$block_size[$i]) . "\t" . $gene_symbol . "\n";
        }
        close OUT;

        last;
      }
    }
  }
}

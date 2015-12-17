#!/usr/bin/perl -w

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);

my ($help, $datafile, $dcolumn, $weightfile, $wcolumn, $minmid, $maxmid, $nbin);
my ($datalist, $weightlist);
my ($tmpx, $i);

$dcolumn = 1;
$wcolumn = 1;

GetOptions('help|?|h' => \$help, 
		   'datafile|f=s' => \$datafile, 
		   'dcolumn|i=i' => \$dcolumn, 
		   'weightfile|w=s' => \$weightfile, 
		   'wcolumn|j=i' => \$wcolumn, 
		   'minmid|n=f' => \$minmid, 
		   'maxmid|x=f' => \$maxmid,
		   'nbin|b=i' => \$nbin);
&usage() if ((! defined($datafile)) or defined($help));

# read data, and decide minimum and maximum bin coordinates of histogram

if (defined($minmid) && defined($maxmid)) {
	($datalist, $tmpx) = &read_column_from_file($datafile, $dcolumn, 0);
}
else {
	($datalist, $tmpx) = &read_column_from_file($datafile, $dcolumn, 1);
	$minmid = $tmpx->[0];
	$maxmid = $tmpx->[1];
}

# read weight
if (defined($weightfile)) {
	($weightlist, $tmpx) = &read_column_from_file($weightfile, $wcolumn, 0);
}
else {
	$weightlist = [(1.0)x($#{$datalist}+1)];
}

# define bin coordinates:
my $binwidth = ($maxmid - $minmid)*1.0/($nbin - 1);
my @binmid = ($minmid);
my @binborders = ();
foreach $i (1..($nbin-1)) {
	$tmpx = $minmid + $i*$binwidth;
	push @binmid, $tmpx;
	$tmpx -= 0.5*$binwidth;
	push @binborders, $tmpx;
}

# calculate the histogram
my @problist = (0.0)x$nbin;
my $norm = 0.0;
foreach $i (0..$#{$datalist}) {
	$tmpx = &find_bin_4_datum($datalist->[$i], \@binborders);
	$problist[$tmpx] += $weightlist->[$i];
	$norm += $weightlist->[$i];
}
# print out
foreach $i (0..($nbin-1)) {
	$problist[$i] /= ($norm*$binwidth);
	printf("%12g\t%12g\n", $binmid[$i], $problist[$i]);
}



####################################################################################

sub usage {
  print "Unknown option: @_\n" if ( @_ );
  print "usage: program --datafile|-f FILENAME [--dcolumn|-i COLUMN] [--weightfile|-w FILENAME] [--wcolumn|-j COLUMN] 
[--minmid|-n MINMID] [--maxmid|-x MAXMID] --nbin|-b NBIN [--help|-?|-h]\n";
  print " --datafile|-f     data file name.\n";
  print " --dcolumn|-i      data column, start from one, not zero.\n";
  print " --weightfile|-w   weight file name.\n";
  print " --wcolumn|-j      wieght column, start form one.\n";
  print " --minmid|-n       the coodinate of the first bin.\n";
  print " --maxmid|-x       the coodinate of the last bin.\n";
  print " --nbin|-b         the number of bins.\n";
  print " --help|-?|-h      this help message.\n";
  exit;
}

sub read_column_from_file {
    my ($filename, $icolumn, $checkmx) = @_; # if $icolumn==0, read the whole line
    local (*INPUT);
	my $line;
	my @split_dummy;
    my @readarray = ();
	my @minmax = ();
    $icolumn -= 1;
    open INPUT, "< $filename";
    while (defined($line=<INPUT>)) {
		unless ($line =~ /^#/) {
			chomp($line);
			$line =~ s/^\s+|\s+$//g; #remove spaces at the beginning or the end of the string
			@split_dummy = split(/\s+/, $line);
			push(@readarray, $split_dummy[$icolumn]);
			if ($checkmx == 1) {
				if ($#readarray == 0) {
					$minmax[0] = $readarray[0];
					$minmax[1] = $readarray[0];
				}
				else {
					$minmax[0] = ($minmax[0], $split_dummy[$icolumn])[$minmax[0] > $split_dummy[$icolumn]];
					$minmax[1] = ($minmax[1], $split_dummy[$icolumn])[$minmax[1] < $split_dummy[$icolumn]];
				}
			} 
		}
	}
	close INPUT;
	return(\@readarray, \@minmax); 
}


sub find_bin_4_datum {
# find out which region a point belongs to 
    my ($value, $borderarray) = @_; # $borderarry is a reference
	my $iregion = -1;

	my $minindex = 0;
	my $maxindex = $#{$borderarray};
	my $midindex;

	if ($value <= $borderarray->[$minindex]) {
		$iregion = 0;
	}
	elsif ($value > $borderarray->[$maxindex]) {
		$iregion = $maxindex + 1;
	}
	else {
		while (($maxindex - $minindex) != 1) {
			$midindex = int(($maxindex + $minindex)/2 + 0.5);
			if ($borderarray->[$midindex] >= $value) {
				$maxindex = $midindex;
			}
			else {
				$minindex = $midindex;
			}
		}
		$iregion = $minindex + 1;
	}
	return($iregion);
}

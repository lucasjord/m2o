#!/usr/bin/perl -w

# Perl script for extracting scans from the DBBC3 system

# Standard perl modules
use Carp;
use Getopt::Long;
use strict;

# Non-standard perl modules
use Astro::Time;

my $pi = $Astro::Time::PI;

# Command line options
my $help=0;
my $input='scan';
my $nbbc=48;
my $output='dbbc_scan.m';
my $version=0;

GetOptions('help!'=>\$help,
	   'input=s'=>\$input,
	   'nbbc=i'=>\$nbbc,
	   'output=s'=>\$output,
           'version!'=>\$version);

# Global variables
my $VERSION=1.1;
my $DATE='12 December 2018';

# If the user requires it print some help
if ($help) {
  print <<HELP;
  This perl script is intended to extract scan data from the DBBC3 system

      Usage : $0 [options]

    The following options are available
      --help                 Produce this output.
      --input=<name>         The root name of the logfiles containing the 
                             dbbc output data. The log files should be
                             <name>_dbbc01.log etc (default=scan)
      --nbbc=<number BBC>    The number of BBCs that there are logfiles for 
                             (default=48).
      --output=<filename>    A MATLAB .m output file for further processing
                             of scan data in matlab (default=dbbc_scan.m)
      --version              Give the version number of the script.
HELP
} elsif ($version) {
  print <<VERSION;
  $0 ; Version : $VERSION ; Date : $DATE
VERSION
}

if ($help || $version) {exit(0);}

my (@azpt1, @elpt1, @azoff1, @eloff1, @power1, @rms1, @offset1, @freqs);
my (@azpt2, @elpt2, @azoff2, @eloff2, @power2, @rms2, @offset2);

# Go through each of the BBC files and extract the scans
for (my $i=1 ; $i <= $nbbc ; $i++) {
  my $file = sprintf "%s_dbbc%02d.log", $input,$i;
  my $azscan = 0;
  open (INPUT, $file) or
    croak "Problem opening scan log file $file ($!)\n";
  my (@fazpt, @felpt, @fazoff, @feloff, @fpower, @frms, @foffset);
  my (@sazpt, @selpt, @sazoff, @seloff, @spower, @srms, @soffset);
  while (<INPUT>) {
    if (/^# Freq : (\S+)\s+MHz/) {
      # Comment line with the frequency information
      push @freqs, $1;
      $azscan = 1;
    } elsif (/^#/) {
      # Comment line, ignore
    } elsif (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)
	      \s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+
	     (\S+)\s+(\S+)\s+(\S+)\s*/x) {
      my $mjd = $1;
      my $ut = $2;
      my $day = $3;
      my $month = $4;
      my $year = $5;
      my $time = $6;
      my $tint = $7;
      if ($azscan) {
	push @fazpt,$8;
	push @felpt,$9;
	push @foffset,$10;
	push @fazoff,$11;
	push @feloff,$12;
	push @fpower,$13;
	push @frms,$14;
      } else {
	push @sazpt,$8;
	push @selpt,$9;
	push @soffset,$10;
	push @sazoff,$11;
	push @seloff,$12;
	push @spower,$13;
	push @srms,$14;
      }
      my $srcname = $15;
      my $stype = $16;
    } elsif (/^\s*$/) {
      if ($azscan) {
	printf "%d points in azimuth scan for file %s\n", scalar @fpower, $file;
      }
      $azscan = 0;
    }
  }
  printf "%d points in elevation scan for file %s\n", scalar @spower, $file;
  close INPUT;
  push @azpt1, \@fazpt;
  push @elpt1, \@felpt;
  push @azoff1, \@fazoff;
  push @eloff1, \@feloff;
  push @power1, \@fpower;
  push @rms1, \@frms;
  push @offset1, \@foffset;
  push @azpt2, \@sazpt;
  push @elpt2, \@selpt;
  push @azoff2, \@sazoff;
  push @eloff2, \@seloff;
  push @power2, \@spower;
  push @rms2, \@srms;
  push @offset2, \@soffset;
}

if (defined $output) {
  open (OUTPUT, ">$output") or
    croak "Problem opening matlab output file $output ($!)\n";
  print OUTPUT "clear all;\n";
  print OUTPUT "close all;\n";
  print OUTPUT "freq = [";
  for (my $i=0 ; $i<scalar @freqs ; $i++) {
    printf OUTPUT "%i ",$freqs[$i];
  }
  print OUTPUT "];\n";      
  my $rows = @offset1;
  my $cols = @{$offset1[0]};
  print OUTPUT "off1 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$offset1[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  print OUTPUT "rms1 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$rms1[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  print OUTPUT "power1 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$power1[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  print OUTPUT "off2 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$offset2[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  print OUTPUT "rms2 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$rms2[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  print OUTPUT "power2 = [";
  for (my $i=0 ; $i<$rows ; $i++) {
    for (my $j=0 ; $j<$cols ; $j++) {
      printf OUTPUT "%.3f ",$power2[$i][$j];
    }
    if ($i < $rows-1) {
      print OUTPUT ";\n";
    } else {
      print OUTPUT "];\n";      
    }
  }
  # Put some code in to make some basic plots.
  print OUTPUT "cut=2.0;\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off1(i,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms1(i,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms1(i,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off1(i,indx),power1(i,indx),'kx');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('Az offset (deg)');\n";
  print OUTPUT "end\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off2(i,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms2(i,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms2(i,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off2(i,indx),power2(i,indx),'ro');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('El offset (deg)');\n";
  print OUTPUT "end\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off1(i+16,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms1(i+16,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms1(i+16,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off1(i+16,indx),power1(i+16,indx),'kx');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i+16));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('Az offset (deg)');\n";
  print OUTPUT "end\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off2(i+16,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms2(i+16,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms2(i+16,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off2(i+16,indx),power2(i+16,indx),'ro');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i+16));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('El offset (deg)');\n";
  print OUTPUT "end\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off1(i+32,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms1(i+32,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms1(i+32,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off1(i+32,indx),power1(i+32,indx),'kx');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i+32));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('Az offset (deg)');\n";
  print OUTPUT "end\n";
  print OUTPUT "figure;\n";
  print OUTPUT "for i = 1:16\n";
  print OUTPUT "npt = length(off2(i+32,:));\n";
  print OUTPUT "indx = true(npt,1);\n";
  print OUTPUT "thres = cut*median(rms2(i+32,:));\n";
  print OUTPUT "for j = 1:npt\n";
  print OUTPUT "if (rms2(i+32,j) > thres)\n";
  print OUTPUT "indx(j) = false;\n";
  print OUTPUT "end\n";
  print OUTPUT "end\n";
  print OUTPUT "subplot(4,4,i);\n";
  print OUTPUT "plot(off2(i+32,indx),power2(i+32,indx),'ro');\n";
  print OUTPUT "lab = sprintf('Freq: %i MHz',freq(i+32));\n";
  print OUTPUT "title(lab);\n";
  print OUTPUT "xlabel('El offset (deg)');\n";
  print OUTPUT "end\n";
  close OUTPUT;
}
